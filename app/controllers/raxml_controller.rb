class RaxmlController < ApplicationController
  
  def index
    @dna_model_options = ""
    models = ["GTRGAMMA","GTRCAT","GTRCAT_FLOAT", "GTRCATI", "GTRGAMMA_FLOAT","GTRGAMMAI"]
    models.each {|m| @dna_model_options= @dna_model_options+"<option>#{m}</option>"}
    @aa_model_options = ""
    models = ["PROTGAMMA","PROTGAMMAI","PROTCAT","PROTCATI"]
    models.each {|m| @aa_model_options= @aa_model_options+"<option>#{m}</option>"}
    @aa_matrices = ""
    matrices = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "GTR"]
    matrices.each {|m| @aa_matrices= @aa_matrices+"<option>#{m}</option>"}
    @heuristics =""
    heuristics = ["MP","ML"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    @heuristics_values =""
    heuristics_values = ["1/2","1/4","1/8","1/16","1/32","1/64"]
    heuristics_values.each {|h| @heuristics_values = @heuristics_values+"<option>#{h}</option>"}
    @raxml = Raxml.new
  end

  def submitJob
    
    @dna_model_options = ""
    models = ["GTRGAMMA","GTRCAT","GTRCAT_FLOAT", "GTRCATI", "GTRGAMMA_FLOAT","GTRGAMMAI"]
    models.each {|m| @dna_model_options= @dna_model_options+"<option>#{m}</option>"}
    @aa_model_options = ""
    models = ["PROTGAMMA","PROTGAMMAI","PROTCAT","PROTCATI"]
    models.each {|m| @aa_model_options= @aa_model_options+"<option>#{m}</option>"}
    @aa_matrices = ""
    matrices = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "GTR"]
    matrices.each {|m| @aa_matrices= @aa_matrices+"<option>#{m}</option>"}
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
    if @query.eql?("DNA")
      @substmodel = "#{params[:dna_substmodel]}"
    else
      @substmodel = params[:aa_substmodel]
      @matrix = params[:matrix]
      @sm_float = params[:sm_float]
      @substmodel = "#{@substmodel}_#{@matrix}#{@sm_float}"
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

    @raxml = Raxml.new({ :alifile =>params[:raxml][:alifile] , :query => @query, :outfile => @outfile, :speed => @speed, :substmodel => @substmodel, :heuristic => @heuristic, :treefile => params[:treefile][:file], :email => @email, :h_value => @h_value, :errorfile => "", :use_heuristic => @use_heuristic, :use_bootstrap => @use_bootstrap, :b_random_seed => @b_random_seed, :b_runs => @b_runs })
    
    
    if @raxml.save
      buildJobDir(@raxml)
      @alifile = saveInfile(@raxml.alifile, "alignment_file")
      @raxml.update_attribute(:alifile,@alifile)
      @raxml.update_attribute(:outfile,"#{@raxml.id}")
#      @raxml.update_attribute(:outfile,@directory+"results.txt")

      @treefile = saveInfile(@raxml.treefile, "tree_file")
      @raxml.update_attribute(:treefile,@treefile)
      
      link = url_for :controller => 'raxml', :action => 'results', :id => @raxml.id
      @raxml.execude(link,@raxml.id.to_s)
      sleep 2
      redirect_to :action => 'wait', :id => @raxml.id 
    else
      @raxml.errors.each do |field, error|
        puts field
        puts error
      end
      render :action => 'index'
    end
  end

  def buildJobDir(rax)
    @directory = "#{RAILS_ROOT}/public/jobs/#{rax.id}/"
    Dir.mkdir(@directory) rescue system("rm -r #{@directory}; mkdir #{@directory}")
  end

  def saveInfile(stream, file_name)
    file_name = @directory+file_name
    #file = @directory+stream.original_filename
    File.open(file_name, "wb") { |f| f.write(stream) }
    return file_name
  end
   
  def wait
    rax = Raxml.find_by_id(params[:id])
   # if pidAlive?(rax.pid)
    if !(jobIsFinished?(params[:id]))
      render :action => "wait"
    else
      redirect_to :action => "results" , :id => rax.id
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
      if f.readlines.size > 1
        if file =~ /submit\.sh\.e/
          @raxml.update_attribute(:errorfile,file)
        end
        finished = true
      end
      f.close
    }
    return finished       
  end

  def results
    rax = Raxml.find_by_id(params[:id])    
    res  =  RaxmlResultsParser.new(rax.outfile)
    @files = res.files
    @names = res.names
    if !(rax.errorfile.eql?(""))
      @files << rax.errorfile
      @names << "error_log"
    end
  end

  def download 
    file = params[:file]
    send_file file
  end
end
