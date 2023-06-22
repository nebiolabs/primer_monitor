class AddDateTriggerToPangolinCalls < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      CREATE FUNCTION init_dates_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.created_at := NOW();
          NEW.updated_at := NOW();
        END;
      $$ LANGUAGE plpgsql;

      CREATE FUNCTION update_date_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.updated_at := NOW();
        END;
      $$ LANGUAGE plpgsql;

      CREATE TRIGGER init_dates BEFORE INSERT ON pangolin_calls
        FOR EACH ROW EXECUTE FUNCTION init_dates_for_pangolin_calls();

      CREATE TRIGGER update_date BEFORE UPDATE ON pangolin_calls
        FOR EACH ROW EXECUTE FUNCTION update_date_for_pangolin_calls();
    SQL
  end

  def down
    execute <<-SQL
      DROP TRIGGER init_dates ON pangolin_calls;
      DROP TRIGGER update_date ON pangolin_calls;

      DROP FUNCTION init_dates_for_pangolin_calls();
      DROP FUNCTION update_date_for_pangolin_calls();
    SQL
  end
end
