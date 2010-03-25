require 'raxml_alignmentfile_parser'
require 'raxml_treefile_parser'
require 'raxml_partitionfile_parser'
require 'raxml_queryfile_parser'


class Raxml < ActiveRecord::Base

  validates_presence_of :alifile, :treefile
  validates_format_of :email, :with => /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i , :on => :save, :message => "Invalid address", :allow_blank => true
  validates_numericality_of :b_random_seed,  :only_integer => true, :greater_than => 0,  :message => "Input must be an integer and greater than 0 and less than 101"
  validates_numericality_of :b_runs, :only_integer => true, :greater_than => 0, :less_than => 101, :message => "Input must be an integer and greater than 0 and less than 101"
  
  validates_presence_of :parfile, :if => :par_selected?
  validates_presence_of :queryfile, :if => :queryfile_selected?
 # validates_each :heuristic  do |record, attr, value|
 #   if (value =~ /^\w+\_([\d+])$/ || value =~ /^\w+\_([\d+]\.{0,1}\d+)$/ )
 #     record.errors.add attr, 'value has to be between 0 and 1' if ($1.to_f > 1.0 || $1.to_f < 0.0)
 #   elsif value =~ /^\w+\_.+$/
 #     record.errors.add attr, 'value is not numeric'
 #   elsif !(value =~ /^none$/) && value =~ /\w+\_$/ 
 #    record.errors.add attr, 'value can\'t be blank when heuristic procedure is selected'
 #   end
 # end
  def par_selected?
    if self.query.eql?("PAR")
      return true
    else
      return false
    end
  end

  def queryfile_selected?
    if self.use_queryfile.eql?("T")
      return true
    else
      return false
    end
  end

  def validate
    if !(self.alifile.eql?("")) && !(self.treefile.eql?(""))
      f = RaxmlAlignmentfileParser.new(self.alifile)
      self.alifile = f.data
      errors.add(:alifile, f.error) if !(f.valid_format)
      t = RaxmlTreefileParser.new(self.treefile)
      self.treefile = t.data
      errors.add(:treefile, t.error) if !(t.valid_format)
      if self.query.eql?("PAR") && !(self.parfile.eql?("") )
        p = RaxmlPartitionfileParser.new(self.parfile,f.ali_length)
        self.parfile = p.data
        errors.add(:parfile, p.error) if !(p.valid_format)
      end
      if self.use_queryfile.eql?("T") && !(self.queryfile.eql?(""))
        q = RaxmlQueryfileParser.new(self.queryfile)
        self.queryfile = q.data
        errors.add(:queryfile, q.error) if !(q.valid_format)
      end
    end
  end

  def execude(link,id)
    #RAxML command with parameters
    
    opts = {"-s" => self.alifile, "-n" => self.outfile, "-m" => self.substmodel,  "-f" => self.speed  , "-link" => link, "-id" => id}
    if emailValid?
      opts["-email"] = self.email
    end
    if !(self.treefile.nil?)
      opts["-t"] = self.treefile
    end
    if self.query.eql?("PAR")
      opts["-q"] = self.parfile
    end
    if self.use_heuristic.eql?("T")
      if self.heuristic.eql?("MP")
        if self.h_value =~ /(1)\/(\d+)/
          opts["-H"] = (($1.to_f)/($2.to_f)).to_s
        end
      elsif self.heuristic.eql?("ML")
        if self.h_value =~ /(1)\/(\d+)/
          opts["-G"] = (($1.to_f)/($2.to_f)).to_s
        end
      end
    elsif self.use_bootstrap.eql?("T")
      opts["-x"] = self.b_random_seed
      opts["-N"] = self.b_runs
    end
    if self.use_queryfile.eql?("T")
      opts["-useQ"] = self.queryfile
    end
    if self.use_clustering.eql?("T")
      opts["-useCl"] = self.use_clustering
    end


    path = "#{RAILS_ROOT}/public/jobs/#{id}"
    shell_file = "#{RAILS_ROOT}/public/jobs/#{id}/submit.sh"
    command = "#{RAILS_ROOT}/bioprogs/ruby/raxml_and_send_email.rb"
    opts.each_key {|k| command  = command+" "+k+" #{opts[k]} "}
    puts command
    File.open(shell_file,'wb'){|file| file.write(command)}
    system "qsub -o #{path} -e #{path} #{shell_file}"
#    process = fork {system command}
#    pid = process+1
#    self.update_attribute(:pid,pid)
  end

  def emailValid?
    if self.email =~ /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i
      return true
    else
      return false
    end      
  end

end
