class CreateUserinfos < ActiveRecord::Migration
  def self.up
    create_table :userinfos do |t|
      t.string :ip,:limit => 16,  :null => false 
      t.integer :saved_submissions, :overall_submissions
      t.timestamps
    end
    add_index(:userinfos, :ip ,:unique => true)
  end

  def self.down
    drop_table :userinfos
  end
end
