class CreateRaxmls < ActiveRecord::Migration
  def self.up
    create_table :raxmls do |t|

      t.string :query, :alifile, :treefile, :outfile, :substmodel, :speed, :heuristic, :h_value, :email ,:errorfile ,:parfile , :queryfile 
      t.integer :b_random_seed, :b_runs 
      t.string :use_heuristic, :use_bootstrap , :use_queryfile, :use_clustering, :default => 'F' , :limit => 1
      t.timestamps
    end
  end

  def self.down
    drop_table :raxmls
  end
end
