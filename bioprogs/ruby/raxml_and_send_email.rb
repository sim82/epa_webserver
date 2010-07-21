#!/usr/bin/ruby

RAILS_ROOT = File.expand_path(File.join(File.dirname(__FILE__), '../..'))
require 'net/smtp'
require "#{RAILS_ROOT}/bioprogs/ruby/reformat.rb"
require "#{File.dirname(__FILE__)}/../../config/environment.rb"
SERVER_NAME = ENV['SERVER_NAME']

class RaxmlAndSendEmail 

  def initialize(opts)
    @raxml_options = Hash.new
    @email_address = ""
    @queryfile = ""
    @use_queryfile = false
    @use_clustering = false
    @link = ""
    @id = ""
    options_parser!(opts)
    @jobpath = "#{RAILS_ROOT}/public/jobs/#{@id}/"
    if @use_clustering
      useUClust
    end
    if @use_queryfile
      buildAlignmentWithHMMER
    end
    exit(0)
    run_raxml
    convertTreefileToPhyloXML
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
      elsif opts[i].eql?("-x")
        @raxml_options["-x"] = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-N")
        @raxml_options["-N"] = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-q")
        @raxml_options["-q"] = opts[i+1]
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
      elsif opts[i].eql?("-useQ")
        @use_queryfile = true
        @queryfile = opts[i+1]
        i = i+1
      elsif opts[i].eql?("-useCl")
        @use_clustering = true
        i = i+1
      else
        raise "ERROR in options_parser!, unknown option #{opts[i]}!"
      end
      i = i+1
    end
  end

  def useUClust
    outfile = @jobpath+"cluster"
    command = "#{RAILS_ROOT}/bioprogs/uclust/uclust32 --input #{@queryfile}  --uc #{outfile}.uc --id 0.90 --usersort 2>&1"
    puts command
    system command
    command = "#{RAILS_ROOT}/bioprogs/uclust/uclust32 --uc2fasta #{outfile}.uc --input #{@queryfile} --output #{outfile}.fas  2>&1"
    puts command
    system command
    ref = Reformat.new(outfile+".fas")
    ref.exportClusterRepresentatives!
    ref.writeToFile(@queryfile)
  end

  def buildAlignmentWithHMMER
    ref = Reformat.new(@raxml_options["-s"])
    ref.reformatToStockholm
    ref.writeToFile(@jobpath+"alignment_file.sto")
    command = "hmmbuild   #{@jobpath}alignment_file.hmm #{@jobpath}alignment_file.sto"
    puts command
    system command
    command = "hmmalign -o #{@jobpath}alignment_file2.sto --mapali #{@jobpath}alignment_file.sto  #{@jobpath}alignment_file.hmm #{@queryfile}  "
    puts command
    system command
    ref = Reformat.new("#{@jobpath}alignment_file2.sto")
    ref.reformatToPhylip
    ref.writeToFile(@raxml_options["-s"])
  end

  def run_raxml
    command ="cd #{RAILS_ROOT}/public/jobs/#{@id}; #{RAILS_ROOT}/bioprogs/raxml/raxmlHPC "
    @raxml_options.each_key  {|k| command = command + k + " " + @raxml_options[k] + " "}
    system command 
  end

  def convertTreefileToPhyloXML
    treefile = @raxml_options["-t"]
    command = "cd #{RAILS_ROOT}/bioprogs/java; java -jar convertToPhyloXML.jar #{@jobpath}RAxML_originalLabelledTree.#{@id}"
    if @raxml_options["-x"].nil? # bootstrapping activated?
      file = @jobpath+"RAxML_classificationLikelihoodWeights."+@id
      command = "#{command} #{file} #{@raxml_options["-N"]} > #{@jobpath}treefile.phyloxml"
      puts command
      system command
      return
    else
      file  = @jobpath+"RAxML_classification."+@id
      command = "#{command} #{file} > #{@jobpath}treefile.phyloxml"
      puts command
      system command
      return
    end
  end

  def send_email
    Net::SMTP.start('localhost', 25) do |smtp|
      smtp.open_message_stream("#{ENV['SERVER_NAME']}", @email_address) do |f|
        
        f.puts "From: #{ENV['SERVER_NAME']}"

        f.puts "To: #{@email_address}"

        f.puts 'Subject: Your RAxML job has been finished.'

        f.puts "Check your results here: #{@link}"

      end

    end
  end
end

RaxmlAndSendEmail.new(ARGV)
