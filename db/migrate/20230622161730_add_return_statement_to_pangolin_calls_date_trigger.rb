class AddReturnStatementToPangolinCallsDateTrigger < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      CREATE OR REPLACE FUNCTION init_dates_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.created_at := NOW();
          NEW.updated_at := NOW();
          RETURN NEW;
        END;
      $$ LANGUAGE plpgsql;

      CREATE OR REPLACE FUNCTION update_date_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.updated_at := NOW();
          RETURN NEW;
        END;
      $$ LANGUAGE plpgsql;
    SQL
  end

  def down
    execute <<-SQL
      CREATE OR REPLACE FUNCTION init_dates_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.created_at := NOW();
          NEW.updated_at := NOW();
        END;
      $$ LANGUAGE plpgsql;

      CREATE OR REPLACE FUNCTION update_date_for_pangolin_calls() RETURNS trigger as $$
        BEGIN
          NEW.updated_at := NOW();
        END;
      $$ LANGUAGE plpgsql;
    SQL
  end
end
