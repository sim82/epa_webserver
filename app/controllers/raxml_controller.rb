class RaxmlController < ApplicationController
  def index
    @dna_model_options = ""
    models = ["GTRGAMMA","GTRCAT", "GTRCATI","GTRGAMMAI"]
    models.each {|m| @dna_model_options= @dna_model_options+"<option>#{m}</option>"}
    @aa_model_options = ""
    models = ["PROTGAMMA","PROTGAMMAI","PROTCAT","PROTCATI"]
    models.each {|m| @aa_model_options= @aa_model_options+"<option>#{m}</option>"}
    @aa_matrices = ""
    matrices = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG"]
    matrices.each {|m| @aa_matrices= @aa_matrices+"<option>#{m}</option>"}
    @par_model_options  =""
    models = ["GAMMA", "GAMMAI", "CAT", "CATI"]
    models.each {|m| @par_model_options= @par_model_options+"<option>#{m}</option>"}
    @heuristics =""
    heuristics = ["MP","ML"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    @heuristics_values =""
    heuristics_values = ["1/2","1/4","1/8","1/16","1/32","1/64"]
    heuristics_values.each {|h| @heuristics_values = @heuristics_values+"<option>#{h}</option>"}
    @raxml = Raxml.new
  end

  def submitJob
    @jobid = generateJobID
    puts "####{@jobid}"
    @dna_model_options = ""
    models = ["GTRGAMMA","GTRCAT", "GTRCATI","GTRGAMMAI"]
    models.each {|m| @dna_model_options= @dna_model_options+"<option>#{m}</option>"}
    @aa_model_options = ""
    models = ["PROTGAMMA","PROTGAMMAI","PROTCAT","PROTCATI"]
    models.each {|m| @aa_model_options= @aa_model_options+"<option>#{m}</option>"}
    @aa_matrices = ""
    matrices = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG"]
    matrices.each {|m| @aa_matrices= @aa_matrices+"<option>#{m}</option>"}
    @par_model_options  =""
    models = ["GAMMA", "GAMMAI", "CAT", "CATI"]
    models.each {|m| @par_model_options= @par_model_options+"<option>#{m}</option>"}
    @heuristics =""
    heuristics = ["MP","ML"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    @heuristics_values =""
    heuristics_values = ["1/2","1/4","1/8","1/16","1/32","1/64"]
    heuristics_values.each {|h| @heuristics_values = @heuristics_values+"<option>#{h}</option>"}
    
    @direcrory = nil
    
    @query = params[:query]
    @speed = params[:speed]
    @substmodel = ""
    @matrix = nil
    @sm_float = nil
    @parfile = ""
    if @query.eql?("DNA")
      @substmodel = "#{params[:dna_substmodel]}"
    elsif @query.eql?("AA")
      @substmodel = params[:aa_substmodel]
      @matrix = params[:matrix]
      @sm_float = params[:sm_float]
      @substmodel = "#{@substmodel}_#{@matrix}#{@sm_float}"
    elsif @query.eql?("PAR")
      @substmodel = "GTR#{params[:par_substmodel]}"
      @parfile = params[:raxml][:parfile]
    end
    @use_clustering = params[:cluster]
    if !(@use_clustering.eql?("T"))
      @use_clustering ="F"
    end
    @use_heuristic = params[:chHeu]
    @heuristic = ""
    @h_value =""
    if  @use_heuristic.eql?("T")
      @heuristic = params[:heuristic]
      @h_value = params[:heu_float]
    else
      @use_heuristic = "F"
    end
    @b_random_seed = 1234
    @b_runs = 100
    @use_bootstrap = params[:chBoot]
    if @use_bootstrap.eql?("T")
      @b_random_seed = params[:random_seed]
      @b_runs = params[:runs]
    else
      @use_bootstrap = "F"
    end
      
    @email = params[:rax_email]
    @outfile = ""
    @alifile = ""
    @queryfile = ""
    @use_queryfile = params[:qfile]
    if @use_queryfile.eql?("T")
      @queryfile = params[:raxml][:queryfile]
    else
      @use_queryfile = "F"
    end
    buildJobDir
    @raxml = Raxml.new({ :alifile =>params[:raxml][:alifile] , :query => @query, :outfile => @outfile, :speed => @speed, :substmodel => @substmodel, :heuristic => @heuristic, :treefile => params[:treefile][:file], :email => @email, :h_value => @h_value, :errorfile => "", :use_heuristic => @use_heuristic, :use_bootstrap => @use_bootstrap, :b_random_seed => @b_random_seed, :b_runs => @b_runs , :parfile => @parfile, :use_queryfile => @use_queryfile, :queryfile => @queryfile, :use_clustering => @use_clustering, :jobid => @jobid})
    
    
    if @raxml.save
 #     @alifile = saveInfile(@raxml.alifile, "alignment_file")
 #     @raxml.update_attribute(:alifile,@alifile)
      @raxml.update_attribute(:outfile,"#{@raxml.jobid}")

#      @treefile = saveInfile(@raxml.treefile, "tree_file")
#      @raxml.update_attribute(:treefile,@treefile)
      
#      if @query.eql?("PAR")
#        @parfile = saveInfile(@raxml.parfile, "partition_file")
#        @raxml.update_attribute(:parfile,@parfile)
#      end

#      if @use_queryfile.eql?("T")
#        @queryfile = saveInfile(@raxml.queryfile, "query_file")
#        @raxml.update_attribute(:queryfile,@queryfile)
#      end
      link = url_for :controller => 'raxml', :action => 'results', :id => @raxml.jobid
      @raxml.execude(link,@raxml.jobid.to_s)
      sleep 2
      redirect_to :action => 'wait', :id => @raxml.jobid 
    else
      @raxml.errors.each do |field, error|
        puts field
        puts error
      end
      render :action => 'index'
    end
  end

  def buildJobDir
    @directory = "#{RAILS_ROOT}/public/jobs/#{@jobid}/"
    Dir.mkdir(@directory) rescue system("rm -r #{@directory}; mkdir #{@directory}")
  end

#  def saveInfile(stream, file_name)
#    file_name = @directory+file_name
#    #file = @directory+stream.original_filename
#    File.open(file_name, "wb") { |f| f.write(stream) }
#    return file_name
#  end

  def generateJobID
    id = "#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}"	
    searching_for_valid_id = true
    while searching_for_valid_id
      r = Raxml.find(:first, :conditions => ["jobid = #{id}"])
      if r.nil?
        return id
      end
      id  = "#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}#{rand(9)}"
    end 
    return id
  end
   
  def wait
    @raxml = Raxml.find(:first, :conditions => ["jobid = #{params[:id]}"])
   # if pidAlive?(rax.pid)
    if !(jobIsFinished?(@raxml.jobid))
      render :action => "wait"
    else
      redirect_to :action => "results" , :id => @raxml.jobid
    end
  end

  def pidAlive?(pid)
    Process.kill(0, Integer(pid))
    return true
  rescue 
    return false
  end

  def jobIsFinished?(id)
    path = "#{RAILS_ROOT}/public/jobs/#{id}/"
    finished = false
    Dir.glob(path+"submit.sh.*"){|file|
      f = File.open(file,'r')
      fi = f.readlines
      if fi.size > 0
        if file =~ /submit\.sh\.e/
          
          @raxml.update_attribute(:errorfile,file)
          f.close
          return true
        else
          fi.each do |line|
            if line =~ /\s+ERROR[\s:]\s*/i
              @raxml.update_attribute(:errorfile,file)
              return true
            elsif line =~ /^done!\s*$/
              return true
            end
          end
        end
      end
      f.close
    }
    return finished       
  end

  def results
    rax =  Raxml.find(:first, :conditions => ["jobid = #{params[:id]}"])
    res  =  RaxmlResultsParser.new(rax.outfile)
    @files = res.files
    @names = res.names
#    @root = "http://lxexelixis1:3000"
    @root  = "#{ENV['SERVER_ADDR']}:3000"
#    @root = "http://lxexelixis1.informatik.tu-muenchen.de:3000"
    @path = "/jobs/#{rax.jobid}/"
    if !(rax.errorfile.eql?(""))
      @files << rax.errorfile
      @names << "logfile"
    end
  end

  def download 
    file = params[:file]
    send_file file
  end
end
