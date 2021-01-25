# frozen_string_literal: true
# This file should contain all the record creation needed to seed the database with its default values.
# The data can then be loaded with the rails db:seed command (or created alongside the database with db:setup).


#password should be changed...
admin = User.create_with(first: 'Admin', last:'User',
                         password: 'admin123',password_confirmation: 'admin123',
                         active: true, approved: true, confirmed: true)
            .find_or_create_by!(email: 'primer-monitor-admin@neb.com')

%w[administrator pi operator participant].each do |role_name|
  Role.find_or_create_by!(name: role_name)
end

admin.roles << Role.find_by(name: 'administrator')
