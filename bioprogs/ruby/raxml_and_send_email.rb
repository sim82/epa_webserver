#!/usr/bin/ruby

RAILS_ROOT = File.expand_path(File.join(File.dirname(__FILE__), '../..'))
require 'net/smtp'
require "#{RAILS_ROOT}/bioprogs/ruby/reformat.rb"
require "#{File.dirname(__FILE__)}/../../config/environment.rb"
require "#{RAILS_ROOT}/bioprogs/ruby/phylip_file_parser.rb"
require "#{RAILS_ROOT}/bioprogs/ruby/fasta_file_parser.rb"
SERVER_NAME = ENV['SERVER_NAME']


### Job processing script 
### Options:
### -s     reference_alignmentfile
### -n     outfile
### -H     MP_value
### -G     ML_value
### -t     treefile
### -f     speed 
### -m     substitution_model
### -x     random_seed
### -N     number_of_bootstrap_samples
### -q     partitionfile
### -email email_address
### -link  link_to_the_results
### -id    jobid
### -useQ  flag if a file with unaligned reads is present
### -useCl flag clustering of the unaligned reads should be performed 
### -mga   flag for using the multigene alignment pipeline

class RaxmlAndSendEmail 

  def initialize(opts)
    @raxml_options = Hash.new
    @email_address = ""
    @queryfile = ""
    @use_queryfile = false
    @use_clustering = false
    @multi_gene_alignment = false
    @link = ""
    @id = ""
    @test_mapping = false
    options_parser!(opts)
    @jobpath = "#{RAILS_ROOT}/public/jobs/#{@id}/"
    if @multi_gene_alignment
      multiGeneAlignment
    else
      if @use_clustering
        useUClust
      end
      if @use_queryfile
        buildAlignmentWithHMMER
      end
      runRAxML
      convertTreefileToPhyloXML
      if @email_address  =~ /\A([^@\s])+@((?:[-a-z0-9]+\.)+[a-z]{2,})\Z/i
        send_email
      end
    end
    puts "done!"
  end

  def options_parser!(opts)
    i = 0
    while i < opts.size
      puts opts[i]
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
      elsif opts[i].eql?("-mga")
        @multi_gene_alignment = true
      elsif opts[i].eql?("-test_mapping")
        @test_mapping = true
      else
        raise "ERROR in options_parser!, unknown option #{opts[i]}!"
      end
      i = i+1
    end
  end

  def multiGeneAlignment
    alifile = @raxml_options["-s"]
    treefile = @raxml_options["-t"]
    partitionfile = @raxml_options["-q"]
    puts partitionfile
    puts alifile
    puts treefile
    puts @queryfile
    outname = @raxml_options["-n"]
    @raxml_options.delete("-q")
    id = @raxml_options["-n"]
    read_gene_mapping = Hash.new
    read_gene_score = Hash.new
    query_sequences = FastaFileParser.new(@queryfile).seqs
    # split the multi gene file based on the partition file
    command ="cd #{RAILS_ROOT}/public/jobs/#{@id}; #{RAILS_ROOT}/bioprogs/raxml/raxmlHPC-SSE3 -s #{alifile} -m GTRGAMMA -p 12345  -n #{id}_mga -fs -t #{treefile} -q #{partitionfile}"
    system command 

    # build fasta db file for each Gene and perform swps3
    gene_no = 0
    Dir.glob(alifile+".GENE.*"){|genefile|
      fasta_db_file = File.new(@jobpath+"GENE#{gene_no}_db.fas", 'w')
      gene = PhylipFileParser.new(genefile)
      gene.disalign!
      gene.seqs.each_key do |k|
        gene.removeGaps!(k)
        if gene.seqs[k].size < 1 # sequences with no letter should be ignored 
          gene.seqs.delete(k)
        else
          fasta_db_file.write(">"+k+"\n"+gene.seqs[k]+"\n")
        end
        
      end
      fasta_db_file.close

      # execute swps3 with every query sequence against the geneX sequences
      query_sequences.each_key do |k|
        query_seq = File.new(@jobpath+"single_query.fas",'w')
        query_seq.write(">#{k}\n#{query_sequences[k]}")
        query_seq.close
        swps3_out = "#{@jobpath}swps3.out"
        command = "#{RAILS_ROOT}/bioprogs/swps3/swps3 #{RAILS_ROOT}/bioprogs/swps3/matrices/blosum62.mat #{query_seq.path} #{fasta_db_file.path} | sort -nr | head -n 1 > #{swps3_out} 2>&1" 
        system command
        # collect query readX x GeneY maximum score
        scorefile = File.open(swps3_out,'r')
        lines = scorefile.readlines
        score = 0
        lines.each do |line|
          if line =~ /^(\d+)/
            score = $1.to_i
            if read_gene_mapping[k].nil?
              read_gene_score[k] = score 
              read_gene_mapping[k] = gene_no
            elsif read_gene_score[k] < score
              read_gene_score[k] = score
              read_gene_mapping[k] = gene_no
            end          
          end
        end
        scorefile.close
      end
      gene_no = gene_no+1
    }
       
    if @test_mapping
      total = 0
      wrong = 0
      read_gene_mapping.each_key do |k|
        if k =~ /GENE_(\d+)/
          real = $1
          assigned = read_gene_mapping[k].to_s
          if !real.eql?(assigned)
            wrong = wrong +1
          end
          total = total+1
        end
      end
      score = 100 - (((wrong.to_i) *100)/(total.to_i))
      puts score.to_s+"% right assigned"
      exit(0)
    end
    # build Gene specific query files
    genes_reads = []
    read_gene_mapping.each_key do |k|
      gene = read_gene_mapping[k]
      if genes_reads[gene].nil?
        genes_reads[gene] = []
        genes_reads[gene] << k
      else
        genes_reads[gene] << k
      end
    end
     i = 0
    
    while i < gene_no
      gene_query_file = "#{@jobpath}queryfile_GENE#{i}.fas"
      out = File.open(gene_query_file,'w')
      genes_reads[i].each do |name|
        out.write(">#{name}\n")
        out.write(query_sequences[name]+"\n")
      end
      # build Alignments with Hmmer and run RAxML
      @raxml_options["-s"] = alifile+".GENE.#{i}"
      @queryfile = gene_query_file
      @raxml_options["-n"] = outname+".GENE#{i}"
      out.close
      buildAlignmentWithHMMER
      puts "RAXML"
      runRAxML
      i = i+1
    end

    # Concat result files
    if @raxml_options["-N"].nil?  # no bootstrap samples, concat RAxML Classification Likelihood weights
      command = "cd #{@jobpath}; cat RAxML_classificationLikelihoodWeights.#{outname}.GENE* >  RAxML_classificationLikelihoodWeights.#{outname}"
    else
      command = "cd #{@jobpath}; cat RAxML_classification.#{outname}.GENE* > RAxML_classification.#{outname}"
    end
    puts command
    system command

    # build phyloXML file
    convertTreefileToPhyloXML
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
    #puts (@raxml_options["-s"])
    ref = Reformat.new(@raxml_options["-s"])
    ref.reformatToStockholm
    ref.writeToFile(@jobpath+"alignment_file.sto")
    command = "hmmbuild --dna  #{@jobpath}alignment_file.hmm #{@jobpath}alignment_file.sto "
    puts command
    system command
    command = "hmmalign -o #{@jobpath}alignment_file2.sto --mapali #{@jobpath}alignment_file.sto  #{@jobpath}alignment_file.hmm #{@queryfile}  "
    puts command
    system command
    ref = Reformat.new("#{@jobpath}alignment_file2.sto")
    ref.reformatToPhylip
    ref.writeToFile(@raxml_options["-s"])
  end

  def runRAxML
    command ="cd #{RAILS_ROOT}/public/jobs/#{@id}; #{RAILS_ROOT}/bioprogs/raxml/raxmlHPC-SSE3 "
    @raxml_options.each_key  {|k| command = command + k + " " + @raxml_options[k] + " "}
    puts command
    system command 
  end

  def convertTreefileToPhyloXML
    treefile = ""
    if @multi_gene_alignment
      treefile = "#{@jobpath}RAxML_originalLabelledTree.#{@id}.GENE0"  
    else
      treefile = "#{@jobpath}RAxML_originalLabelledTree.#{@id}"
    end
    command = "cd #{RAILS_ROOT}/bioprogs/java; java -jar convertToPhyloXML.jar #{treefile}"
    if @raxml_options["-x"].nil? # bootstrapping activated?
      file = @jobpath+"RAxML_classificationLikelihoodWeights."+@id
      command = "#{command} #{file}  > #{@jobpath}treefile.phyloxml"
      puts command
      system command
      return
    else
      file  = @jobpath+"RAxML_classification."+@id
      command = "#{command} #{file} #{@raxml_options["-N"]}> #{@jobpath}treefile.phyloxml"
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
