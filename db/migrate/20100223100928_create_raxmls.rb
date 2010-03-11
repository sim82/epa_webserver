class CreateRaxmls < ActiveRecord::Migration
  def self.up
    create_table :raxmls do |t|

      t.string :query, :alifile, :treefile, :outfile, :substmodel, :speed, :heuristic, :h_value, :email 
      t.integer :pid
            
      t.timestamps
    end
  end

  def self.down
    drop_table :raxmls
  end
end
