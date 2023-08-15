# frozen_string_literal: true

class PrimerSet < ApplicationRecord
  belongs_to :user
  belongs_to :organism
  belongs_to :amplification_method
  has_many :oligos, -> { order('oligos.ref_start') }, dependent: :destroy, inverse_of: :primer_set
  has_many :subscriptions, dependent: :destroy, class_name: 'PrimerSetSubscription'
  has_many :subscribers, through: :subscriptions, source: :user

  enum status: { pending: 'pending', complete: 'complete', failed: 'failed' }

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

  after_save :notify_admins_about_primer_set_update

  def to_s
    name
  end

  def display_url
    citation_url || doi_url
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
      PrimerSetMailer.updated_primer_set_email(user.email, self).deliver_later
    end
  end

  def align_primers
    pid = Process.spawn({ 'DB_HOST' => ENV['DB_HOST'], 'DB_NAME' => ENV['DB_NAME'], 'DB_USER' => ENV['DB_USER'] },
                        Shellwords.join(['bash', 'lib/update_primers.sh', 'lib/cov_index/NC_045512.2',
                                         id.to_s]))
    Process.detach pid # prevent zombie process
    pid
  end
end
