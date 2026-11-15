# frozen_string_literal: true

require 'fileutils'

class PrimerSet < ApplicationRecord
  belongs_to :user
  belongs_to :organism
  belongs_to :amplification_method
  has_many :oligos, -> { order('oligos.name') }, dependent: :destroy, inverse_of: :primer_set
  has_many :subscriptions, dependent: :destroy, class_name: 'PrimerSetSubscription'
  has_many :subscribers, through: :subscriptions, source: :user

  enum :status, { created: 'created', complete: 'complete', failed: 'failed', processing: 'processing' }

  # The review lifecycle. A new set waits for an administrator. Approving it sets processing, and the backend cron
  # job (lib/backend_scripts/process-primer-sets.sh) then builds the public tracks for every processing set whose
  # oligos are all aligned and marks it complete. Nothing sets failed automatically; it means rejected.
  STATUS_LABELS = {
    'created' => 'Awaiting review',
    'processing' => 'Approved, publishing',
    'complete' => 'Published',
    'failed' => 'Rejected'
  }.freeze

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

  validates :name, uniqueness: true, presence: true
  HTTP_URL = %r{\Ahttps?://\S+\z}i
  # only when changed, so sets saved before this rule stay editable; display_url doesn't link their URLs
  validates :citation_url, format: { with: HTTP_URL, message: 'must start with http:// or https://' },
                           allow_blank: true, if: :citation_url_changed?

  validates :oligos, presence: true

  validates_associated :oligos

  # after commit, not after save: the :async job adapter runs the mailer job on another thread right away, and
  # before the commit that thread can't see the record yet (ActiveJob::DeserializationError, no email)
  after_save_commit :notify_admins_about_primer_set_update

  def to_s
    name
  end

  def display_url
    (citation_url if citation_url&.match?(HTTP_URL)) || doi_url
  end

  def doi_url
    "https://doi.org/#{doi}" if doi.present?
  end

  def status_label
    STATUS_LABELS.fetch(status.to_s, status.to_s)
  end

  # the backend only publishes a set once every oligo has an alignment position
  def unaligned_oligos_count
    oligos.where.not(id: OligoAlignmentPosition.select(:oligo_id)).count
  end

  def subscription_for_user(user)
    return unless user

    subscriptions.where(user_id: user.id).first
  end

  def notify_admins_about_primer_set_update
    Role.find_by(name: 'administrator').users.where(active: true).find_each do |user|
      PrimerSetMailer.updated_primer_set_email(user.email, self).deliver_later unless user.email == ENV['ADMIN_EMAIL']
    end
  end

  STRANDS = { 'Forward' => '+', 'Reverse' => '-' }.freeze

  # Where the alignment step (lib/update_primers.sh) placed each oligo on the reference, as igv.js features.
  # This lets the show page draw the set before it is approved and published to the data server.
  def alignment_features
    positions = OligoAlignmentPosition.joins(:oligo, :organism_taxon)
                                      .where(oligos: { primer_set_id: id }).order(:ref_start)
    rows = positions.pluck('organism_taxa.reference_accession', :ref_start, :ref_end, 'oligos.name', 'oligos.category')
    rows.map do |chr, start, stop, name, category|
      { chr:, start:, end: stop, name:, strand: STRANDS[category], category: }.compact
    end
  end

  # TODO: switch this to use delayed job, avoid multiple alignments in succession
  def align_primers
    shared_dir = ENV.fetch('DEPLOY_SHARED_DIR', nil)
    index_name = organism.name.parameterize
    # the script also reads DB_HOST, DB_NAME, DB_USER and PIXI_BIN_PATH, which it inherits from this process
    pid = Process.spawn({ 'PGPASSFILE' => "#{shared_dir}/config/.pgpass" },
                        'bash', 'lib/update_primers.sh', 'lib/alignment/pixi.toml',
                        "bt2_indices/#{index_name}/#{index_name}", id.to_s,
                        out: [primer_alignment_log_path, 'a'], err: %i[child out])
    Process.detach pid # prevent zombie process
    pid
  end

  private

  def primer_alignment_log_path
    base_dir = ENV['FRONTEND_LOG_PATH'].presence || Rails.root.join('log').to_s
    FileUtils.mkdir_p(base_dir)
    File.join(base_dir, 'primer_alignment.log')
  rescue SystemCallError
    fallback_dir = Rails.root.join('log').to_s
    FileUtils.mkdir_p(fallback_dir)
    File.join(fallback_dir, 'primer_alignment.log')
  end
end
