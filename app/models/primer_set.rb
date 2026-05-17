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

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

  validates :name, uniqueness: true, presence: true

  validates :oligos, presence: true

  validates_associated :oligos

  after_save :notify_admins_about_primer_set_update

  def to_s
    name
  end

  def display_url
    citation_url.presence || doi_url
  end

  def doi_url
    "https://doi.org/#{doi}" if doi.present?
  end

  def subscription_for_user(user)
    return unless user

    subscriptions.where(user_id: user.id).first
  end

  def notify_admins_about_primer_set_update
    Role.find_by(name: 'administrator').users.each do |user|
      PrimerSetMailer.updated_primer_set_email(user.email, self).deliver_later unless user.email == ENV['ADMIN_EMAIL']
    end
  end

  # TODO: switch this to use delayed job, avoid multiple alignments in succession
  def align_primers
    log_path = primer_alignment_log_path
    pid = Process.spawn({ 'DB_HOST' => ENV['DB_HOST'], 'DB_NAME' => ENV['DB_NAME'], 'DB_USER' => ENV['DB_USER'],
                          'MICROMAMBA_BIN_PATH' => ENV['MICROMAMBA_BIN_PATH'],
                          'PGPASSFILE' => "#{ENV['DEPLOY_SHARED_DIR']}/config/.pgpass" },
                        Shellwords.join(['bash', 'lib/update_primers.sh', "#{ENV['DEPLOY_SHARED_DIR']}/alignment_env",
                                         "bt2_indices/#{organism.name.parameterize}/#{organism.name.parameterize}",
                                         id.to_s]) +
                          " >> #{Shellwords.escape(log_path)} 2>&1")
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
