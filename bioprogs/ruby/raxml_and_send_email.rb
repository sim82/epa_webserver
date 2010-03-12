#!/usr/bin/ruby

RAILS_ROOT = File.expand_path(File.join(File.dirname(__FILE__), '../..'))
require 'net/smtp'

class RaxmlAndSendEmail 

  def initialize(opts)
    @raxml_options = Hash.new
    @email_address = ""
    @link = ""
    @id = ""
    options_parser!(opts)
    run_raxml
    if @email_address  =~ /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i
      send_email
    end
    puts "done!"
  end

  def options_parser!(opts)
    i = 0
    while i < opts.size
      if opts[i].eql?("-s")
        @raxml_options["-s"] = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-n")
        @raxml_options["-n"] = opts[i+1]
        i = i+1  
      elsif opts[i].eql?("-H")
        @raxml_options["-H"] = opts[i+1]
        i = i+1  
      elsif opts[i].eql?("-G")
        @raxml_options["-G"] = opts[i+1]
        i = i+1      
      elsif opts[i].eql?("-t")
        @raxml_options["-t"] = opts[i+1]
        i = i+1  
      elsif opts[i].eql?("-f")
        @raxml_options["-f"] = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-m")
        a = opts[i+1].split("_")
        @raxml_options["-m"] = a.join(" ")
        i = i+1  
      elsif opts[i].eql?("-email")
        @email_address = opts[i+1]
        i = i+1 
      elsif opts[i].eql?("-link")
        @link = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-id")
        @id = opts[i+1]
        i = i+1
      else
        raise "ERROR in options_parser!, unknown option #{opts[i]}!"
      end
      i = i+1
    end
  end

  def run_raxml
    command ="cd #{RAILS_ROOT}/public/jobs/#{@id}; #{RAILS_ROOT}/bioprogs/raxml/raxmlHPC "
    @raxml_options.each_key  {|k| command = command + k + " " + @raxml_options[k] + " "}
#    puts "********************************************************"
 #   puts command
 #   puts "********************************************************"
    system command 
  end

  def send_email
    Net::SMTP.start('localhost', 25) do |smtp|
      smtp.open_message_stream('raxml@lxexelixis1.informatik.tu-muenchen.de', @email_address) do |f|
        
        f.puts 'From: raxml@lxexelixis1.informatik.tu-muenchen.de'

        f.puts "To: #{@email_address}"

        f.puts 'Subject: Your RAxML job has been finished.'

        f.puts "Check your results here: #{@link}"

      end

    end
  end
end

RaxmlAndSendEmail.new(ARGV)
