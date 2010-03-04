class Raxml < ActiveRecord::Base
  
  validates_presence_of :alifile
  validates_format_of :email, :with => /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i , :on => :save, :message => "Invalid email adress", :allow_blank => true
  validates_each :heuristics :allow_blank => 'true' do |record, attr, value|
    if (value =~ /^\w+\s+([01])$/  ||   value =~ /^\w+\s+([01]\.{0,1}\d+)$/)
      record.errors.add attr, 'Heuristic value has to be between 0 and 1' if ($1.to_f > 1.0 || $1.to_f < 0.0)
    else
      record.errors.add attr, 'Heuristic value is not numeric'
    end
  end
  



  def execude(link)
    #RAxML command with parameters
    
    opts = {"-wait" => self.wait.to_s, "-s" => self.alifile, "-n" => self.outfile, "-m" => self.substmodel,  "-f" => self.speed  , "-link" => link}
    if emailValid?
      opts["-email"] = self.email
    end
    if !(self.treefile.nil?)
      opts["-t"] = self.treefile
    end
    if !(self.heuristic.eql?("none"))
      if self.heuristic =~ /(\w{2})\s([0-9\.]+)/
        if $1.eql?("MP")
          opts["-H"] = $2
        elsif $2.eql?("MC")
          opts["-G"] = $2
        end
      end
    end

    command = "#{RAILS_ROOT}/bioprogs/ruby/raxml_and_send_email.rb"
    opts.each_key {|k| command  = command+" "+k+" "+opts[k]+" "}
    puts command
    process = fork {system command}
    pid = process+1
    self.update_attribute(:pid,pid)
  end

  def emailValid?
    if self.email =~ /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i
      return true
    else
      return false
    end      
  end


  def Raxml.isNumeric?(s)
    if Float(s) rescue  return false # not numeric
 
    else
    return true # numeric
  end


end
