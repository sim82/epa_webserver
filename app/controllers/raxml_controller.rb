class RaxmlController < ApplicationController
  def index
    
    @dna_model_options = ""
    @aa_model_options = ""
    @aa_matrices = ""
    @par_model_options  =""
    @heuristics =""
    @heuristics_values =""
    @ip_counter = 0;
    @submission_counter = 0;
    initialize_options
    @raxml = Raxml.new
  end

  def initialize_options
    models = ["GTRGAMMA","GTRCAT", "GTRCATI","GTRGAMMAI"]
    models.each {|m| @dna_model_options= @dna_model_options+"<option>#{m}</option>"}
    models = ["PROTGAMMA","PROTGAMMAI","PROTCAT","PROTCATI"]
    models.each {|m| @aa_model_options= @aa_model_options+"<option>#{m}</option>"}
    matrices = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG"]
    matrices.each {|m| @aa_matrices= @aa_matrices+"<option>#{m}</option>"}
    models = ["GAMMA", "GAMMAI", "CAT", "CATI"]
    models.each {|m| @par_model_options= @par_model_options+"<option>#{m}</option>"}
    heuristics = ["MP","ML"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    heuristics_values = ["1/2","1/4","1/8","1/16","1/32","1/64"]
    heuristics_values.each {|h| @heuristics_values = @heuristics_values+"<option>#{h}</option>"}
    
    ips = Userinfo.find(:all)
    @ip_counter = (ips.size) - 1   # - c.c.c.c
    userinfo  = Userinfo.find(:first, :conditions => {:ip => "c.c.c.c"})
    @submission_counter = userinfo.overall_submissions
    
  end

  def submitJob
    @jobid = generateJobID
    @dna_model_options = ""
    @aa_model_options = ""
    @aa_matrices = ""
    @par_model_options  =""
    @heuristics =""
    @heuristics_values =""
    @ip_counter = 0;
    @submission_counter = 0;
    initialize_options
    
    @direcrory = nil
    @ip = request.env['REMOTE_ADDR']
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
    @raxml = Raxml.new({ :alifile =>params[:raxml][:alifile] , :query => @query, :outfile => @outfile, :speed => @speed, :substmodel => @substmodel, :heuristic => @heuristic, :treefile => params[:treefile][:file], :email => @email, :h_value => @h_value, :errorfile => "", :use_heuristic => @use_heuristic, :use_bootstrap => @use_bootstrap, :b_random_seed => @b_random_seed, :b_runs => @b_runs , :parfile => @parfile, :use_queryfile => @use_queryfile, :queryfile => @queryfile, :use_clustering => @use_clustering, :jobid => @jobid, :user_ip => @ip})
    
    
    if @raxml.save
      @raxml.update_attribute(:outfile,"#{@raxml.jobid}")
      link = url_for :controller => 'raxml', :action => 'results', :id => @raxml.jobid
      @raxml.execude(link,@raxml.jobid.to_s)

      ## save userinfos
      ip = request.env['REMOTE_ADDR']
      if ip.eql?("") || ip.nil?
        ip = "xxx.xxx.xxx.xxx"
      end
      if Userinfo.exists?(:ip => ip)
        userinfo = Userinfo.find(:first, :conditions => {:ip => ip})
        userinfo.update_attribute(:saved_submissions, userinfo.saved_submissions+1)
        userinfo.update_attribute(:overall_submissions, userinfo.overall_submissions+1)
      else
        userinfo = Userinfo.new({:ip => ip, :saved_submissions => 1, :overall_submissions => 1})
        userinfo.save
      end
      #### main counter c.c.c.c   ## faster but accurate?
      counter_ip = "c.c.c.c"
      if Userinfo.exists?(:ip => counter_ip)
        userinfo = Userinfo.find(:first, :conditions => {:ip => counter_ip })
        userinfo.update_attribute(:saved_submissions, userinfo.saved_submissions+1)
        userinfo.update_attribute(:overall_submissions, userinfo.overall_submissions+1)
      else
        userinfo = Userinfo.new({:ip => counter_ip, :saved_submissions => 1, :overall_submissions => 1})
        userinfo.save
      end

      sleep 2 #Without this, an error occurs. Somehow the writing in the database is not fast enough
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
    @ip = request.env['REMOTE_ADDR']
    if !(jobIsFinished?(@raxml.jobid))
      render :action => "wait"
    else
      redirect_to :action => "results" , :id => @raxml.jobid
    end
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
    @root  = "#{ENV['SERVER_ADDR']}:3000"
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
