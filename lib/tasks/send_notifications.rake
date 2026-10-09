# frozen_string_literal: true

namespace :notifications do
  desc 'Sends notifications about primer overlaps'
  task send: :environment do
    new_pns = ProposedNotification.new_proposed_notifications
    Rails.logger.info("Found #{new_pns.size} new proposed notifications")
    new_pns.each(&:save!)
    VerifiedNotification.find_or_create_verified_notifications!.each(&:deliver!)
  end

  desc 'Records current primer overlaps as already notified, emailing only the given addresses on the next send'
  task :baseline, [:email] => :environment do |_, args|
    new_pns = ProposedNotification.new_proposed_notifications
    new_pns.each(&:save!)
    skipped = VerifiedNotification.skip_unsent!(except_emails: [args[:email], *args.extras].compact)
    puts "Recorded #{new_pns.size} new proposed notifications; skipped emailing #{skipped.size} users"
  end
end
