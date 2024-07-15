# This file is auto-generated from the current state of the database. Instead
# of editing this file, please use the migrations feature of Active Record to
# incrementally modify your database, and then regenerate this schema definition.
#
# This file is the source Rails uses to define your schema when running `bin/rails
# db:schema:load`. When creating a new database, `bin/rails db:schema:load` tends to
# be faster and is potentially less error prone than running all of your
# migrations from scratch. Old migrations may fail to apply correctly if those
# migrations use external dependencies or application code.
#
# It's strongly recommended that you check this file into your version control system.

ActiveRecord::Schema[7.0].define(version: 2024_07_15_163630) do
  # These are extensions that must be enabled in order to support this database
  enable_extension "plpgsql"

  # Custom types defined in this database.
  # Note that some types may not work with other database engines. Be careful if changing database.
  create_enum "oligo_category", ["BIP", "FIP", "LB", "LF", "B3", "F3", "Reverse", "Forward", "Probe"]
  create_enum "primer_set_status", ["created", "complete", "failed", "processing"]

  create_table "amplification_methods", force: :cascade do |t|
    t.string "name"
    t.string "description_url"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "blast_hits", force: :cascade do |t|
    t.bigint "oligo_id"
    t.bigint "organism_id"
    t.integer "num_identities"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["oligo_id"], name: "index_blast_hits_on_oligo_id"
    t.index ["organism_id"], name: "index_blast_hits_on_organism_id"
  end

  create_table "detailed_geo_location_aliases", force: :cascade do |t|
    t.string "world", null: false
    t.string "region"
    t.string "subregion"
    t.string "division"
    t.string "subdivision"
    t.string "locality"
    t.string "sublocality"
    t.float "latitude"
    t.float "longitude"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["region", "subregion", "division", "subdivision", "locality", "sublocality", "latitude", "longitude"], name: "alias_full_record", unique: true
  end

  create_table "detailed_geo_locations", force: :cascade do |t|
    t.string "world", null: false
    t.string "region"
    t.string "subregion"
    t.string "division"
    t.string "subdivision"
    t.string "locality"
    t.string "sublocality"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.bigint "detailed_geo_location_alias_id", null: false
    t.index ["detailed_geo_location_alias_id"], name: "index_detailed_geo_locations_on_detailed_geo_location_alias_id"
    t.index ["region", "subregion", "division", "subdivision", "locality", "sublocality"], name: "full_record", unique: true
  end

  create_table "fasta_records", force: :cascade do |t|
    t.string "strain", null: false
    t.string "genbank_accession"
    t.string "gisaid_epi_isl"
    t.date "date_collected"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.string "variant_name"
    t.bigint "detailed_geo_location_id"
    t.date "date_submitted"
    t.text "pangolin_lineage"
    t.text "pangolin_version"
    t.bigint "lineage_call_id"
    t.bigint "pending_lineage_call_id"
    t.bigint "organism_taxon_id"
    t.index ["created_at"], name: "index_fasta_records_on_created_at"
    t.index ["detailed_geo_location_id"], name: "index_fasta_records_on_detailed_geo_location_id"
    t.index ["lineage_call_id"], name: "index_fasta_records_on_lineage_call_id"
    t.index ["organism_taxon_id"], name: "index_fasta_records_on_organism_taxon_id"
    t.index ["pending_lineage_call_id"], name: "index_fasta_records_on_pending_lineage_call_id"
  end

  create_table "genomic_features", force: :cascade do |t|
    t.string "name"
    t.string "feature_type"
    t.integer "ref_start"
    t.integer "ref_end"
    t.bigint "organism_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["organism_id"], name: "index_genomic_features_on_organism_id"
  end

  create_table "lineage_callers", force: :cascade do |t|
    t.string "name", null: false
    t.string "version_specifiers"
    t.string "script_name", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "pending_version_specifiers"
    t.string "dataset_versions"
    t.string "pending_dataset_versions"
  end

  create_table "lineage_calls", force: :cascade do |t|
    t.string "taxon", null: false
    t.bigint "lineage_id"
    t.bigint "lineage_caller_id"
    t.string "metadata"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["lineage_caller_id"], name: "index_lineage_calls_on_lineage_caller_id"
    t.index ["lineage_id"], name: "index_lineage_calls_on_lineage_id"
  end

  create_table "lineages", force: :cascade do |t|
    t.string "name", null: false
    t.string "caller_name"
    t.bigint "organism_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "external_link"
    t.index ["name"], name: "index_lineages_on_name", unique: true
    t.index ["organism_id"], name: "index_lineages_on_organism_id"
  end

  create_table "oligo_alignment_positions", force: :cascade do |t|
    t.bigint "organism_taxon_id", null: false
    t.bigint "oligo_id", null: false
    t.integer "ref_start", null: false
    t.integer "ref_end", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["oligo_id"], name: "index_oligo_alignment_positions_on_oligo_id"
    t.index ["organism_taxon_id"], name: "index_oligo_alignment_positions_on_organism_taxon_id"
  end

  create_table "oligos", force: :cascade do |t|
    t.string "name"
    t.string "sequence"
    t.bigint "primer_set_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "locus"
    t.enum "category", enum_type: "oligo_category"
    t.string "short_name"
    t.string "strand"
    t.index ["primer_set_id"], name: "index_oligos_on_primer_set_id"
  end

  create_table "oligos_test", id: false, force: :cascade do |t|
    t.integer "id"
    t.integer "ref_start"
    t.integer "ref_end"
    t.integer "primer_set_id"
  end

  create_table "organism_taxa", force: :cascade do |t|
    t.string "name", null: false
    t.string "reference_accession", null: false
    t.bigint "organism_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.integer "ncbi_taxon_id"
    t.bigint "lineage_caller_id"
    t.index ["lineage_caller_id"], name: "index_organism_taxa_on_lineage_caller_id"
    t.index ["organism_id"], name: "index_organism_taxa_on_organism_id"
  end

  create_table "organisms", force: :cascade do |t|
    t.string "name", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "alias"
    t.string "slug"
    t.boolean "public"
    t.integer "variant_bed_lookback_days", default: 180, null: false
    t.index ["slug"], name: "index_organisms_on_slug", unique: true
  end

  create_table "primer_set_subscriptions", force: :cascade do |t|
    t.bigint "primer_set_id", null: false
    t.bigint "user_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.boolean "active", default: true, null: false
    t.index ["primer_set_id"], name: "index_primer_set_subscriptions_on_primer_set_id"
    t.index ["user_id", "primer_set_id"], name: "index_primer_set_subscriptions_on_user_id_and_primer_set_id", unique: true
    t.index ["user_id"], name: "index_primer_set_subscriptions_on_user_id"
  end

  create_table "primer_sets", force: :cascade do |t|
    t.string "name"
    t.bigint "user_id"
    t.integer "organism_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.enum "status", default: "created", enum_type: "primer_set_status"
    t.string "citation_url"
    t.string "doi"
    t.bigint "amplification_method_id"
    t.index ["name"], name: "index_primer_sets_on_name", unique: true
    t.index ["organism_id"], name: "index_primer_sets_on_organism_id"
    t.index ["user_id"], name: "index_primer_sets_on_user_id"
  end

  create_table "proposed_notifications", force: :cascade do |t|
    t.bigint "primer_set_id", null: false
    t.bigint "oligo_id", null: false
    t.bigint "verified_notification_id"
    t.integer "coordinate", null: false
    t.float "fraction_variant", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.bigint "subscribed_geo_location_id", null: false
    t.bigint "primer_set_subscription_id", null: false
    t.bigint "user_id", null: false
    t.bigint "detailed_geo_location_alias_id"
    t.index ["detailed_geo_location_alias_id"], name: "index_proposed_notifications_on_detailed_geo_location_alias_id"
    t.index ["oligo_id"], name: "index_proposed_notifications_on_oligo_id"
    t.index ["primer_set_id"], name: "index_proposed_notifications_on_primer_set_id"
    t.index ["primer_set_subscription_id"], name: "index_proposed_notifications_on_primer_set_subscription_id"
    t.index ["subscribed_geo_location_id"], name: "index_proposed_notifications_on_subscribed_geo_location_id"
    t.index ["user_id"], name: "index_proposed_notifications_on_user_id"
    t.index ["verified_notification_id"], name: "index_proposed_notifications_on_verified_notification_id"
  end

  create_table "roles", force: :cascade do |t|
    t.string "name"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["name"], name: "index_roles_on_name", unique: true
  end

  create_table "subscribed_geo_locations", force: :cascade do |t|
    t.bigint "user_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.bigint "detailed_geo_location_alias_id", null: false
    t.index ["detailed_geo_location_alias_id"], name: "tmp"
    t.index ["user_id", "detailed_geo_location_alias_id"], name: "index_subscribed_geo_locs_on_user_and_detailed_geo_loc_alias_id", unique: true
    t.index ["user_id"], name: "index_subscribed_geo_locations_on_user_id"
  end

  create_table "user_roles", force: :cascade do |t|
    t.bigint "user_id", null: false
    t.bigint "role_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["role_id"], name: "index_user_roles_on_role_id"
    t.index ["user_id", "role_id"], name: "index_user_roles_on_user_id_and_role_id", unique: true
    t.index ["user_id"], name: "index_user_roles_on_user_id"
  end

  create_table "users", force: :cascade do |t|
    t.string "first", null: false
    t.string "last", null: false
    t.string "email", null: false
    t.boolean "activated", default: false, null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "login"
    t.string "encrypted_password"
    t.integer "sign_in_count", default: 0, null: false
    t.integer "failed_login_count", default: 0, null: false
    t.datetime "last_request_at", precision: nil
    t.datetime "current_sign_in_at", precision: nil
    t.datetime "last_sign_in_at", precision: nil
    t.string "current_sign_in_ip"
    t.string "last_sign_in_ip"
    t.boolean "active", default: false
    t.boolean "approved", default: false
    t.boolean "confirmed", default: false
    t.boolean "send_primer_updates", default: false, null: false
    t.string "reset_password_token"
    t.datetime "reset_password_sent_at", precision: nil
    t.string "remember_token"
    t.datetime "remember_created_at", precision: nil
    t.string "authentication_token"
    t.string "confirmation_token", limit: 255
    t.datetime "confirmed_at", precision: nil
    t.datetime "confirmation_sent_at", precision: nil
    t.string "unconfirmed_email"
    t.integer "failed_attempts", default: 0, null: false
    t.string "unlock_token"
    t.datetime "locked_at", precision: nil
    t.string "provider"
    t.string "uid"
    t.integer "lookback_days", default: 30, null: false
    t.float "variant_fraction_threshold", default: 0.1, null: false
    t.index ["confirmation_token"], name: "index_users_on_confirmation_token", unique: true
    t.index ["email"], name: "index_users_on_email", unique: true
    t.index ["reset_password_token"], name: "index_users_on_reset_password_token", unique: true
    t.index ["unconfirmed_email"], name: "index_users_on_unconfirmed_email", unique: true
    t.index ["unlock_token"], name: "index_users_on_unlock_token", unique: true
  end

  create_table "variant_sites", force: :cascade do |t|
    t.integer "ref_start", null: false
    t.string "variant_type", null: false
    t.string "variant", null: false
    t.bigint "fasta_record_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.integer "ref_end", null: false
    t.virtual "usable_del_or_snp", type: :boolean, as: "((((variant_type)::text = 'D'::text) OR ((variant_type)::text = 'X'::text)) AND ((variant)::text !~~ '%N%'::text))", stored: true
    t.virtual "usable_insertion", type: :boolean, as: "(((variant_type)::text = 'I'::text) AND ((variant)::text !~~ '%N%'::text))", stored: true
    t.bigint "organism_taxon_id"
    t.index ["fasta_record_id"], name: "index_variant_sites_on_fasta_record_id"
    t.index ["organism_taxon_id"], name: "index_variant_sites_on_organism_taxon_id"
    t.index ["ref_start", "fasta_record_id", "variant_type"], name: "ensure_variant_unique", unique: true
    t.index ["usable_del_or_snp"], name: "variant_sites_usable_del_or_snp_idx"
    t.index ["usable_insertion"], name: "variant_sites_usable_insertion_idx"
  end

  create_table "verified_notifications", force: :cascade do |t|
    t.bigint "user_id", null: false
    t.date "date_sent"
    t.string "status", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["user_id"], name: "index_verified_notifications_on_user_id"
  end

  add_foreign_key "blast_hits", "oligos"
  add_foreign_key "blast_hits", "organisms"
  add_foreign_key "detailed_geo_locations", "detailed_geo_location_aliases"
  add_foreign_key "fasta_records", "detailed_geo_locations"
  add_foreign_key "fasta_records", "lineage_calls"
  add_foreign_key "fasta_records", "lineage_calls", column: "pending_lineage_call_id"
  add_foreign_key "fasta_records", "organism_taxa"
  add_foreign_key "genomic_features", "organisms"
  add_foreign_key "lineage_calls", "lineage_callers"
  add_foreign_key "lineage_calls", "lineages"
  add_foreign_key "lineages", "organisms"
  add_foreign_key "oligo_alignment_positions", "oligos"
  add_foreign_key "oligo_alignment_positions", "organism_taxa"
  add_foreign_key "oligos", "primer_sets"
  add_foreign_key "organism_taxa", "lineage_callers"
  add_foreign_key "organism_taxa", "organisms"
  add_foreign_key "primer_set_subscriptions", "primer_sets"
  add_foreign_key "primer_set_subscriptions", "users"
  add_foreign_key "primer_sets", "amplification_methods"
  add_foreign_key "primer_sets", "organisms"
  add_foreign_key "primer_sets", "users"
  add_foreign_key "proposed_notifications", "detailed_geo_location_aliases"
  add_foreign_key "proposed_notifications", "oligos"
  add_foreign_key "proposed_notifications", "primer_set_subscriptions"
  add_foreign_key "proposed_notifications", "primer_sets"
  add_foreign_key "proposed_notifications", "subscribed_geo_locations"
  add_foreign_key "proposed_notifications", "users"
  add_foreign_key "proposed_notifications", "verified_notifications"
  add_foreign_key "subscribed_geo_locations", "detailed_geo_location_aliases"
  add_foreign_key "subscribed_geo_locations", "users"
  add_foreign_key "user_roles", "roles"
  add_foreign_key "user_roles", "users"
  add_foreign_key "variant_sites", "fasta_records"
  add_foreign_key "variant_sites", "organism_taxa"
  add_foreign_key "verified_notifications", "users"
end
