class RaxmlController < ApplicationController
  def submit
    @root  = "#{ENV['SERVER_ADDR']}:3000"
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
    
    getInfo
    
  end

  def getInfo
    ips = Userinfo.find(:all)
    if ips.size == 0
      @ip_counter=0;
       @submission_counter = 0;
    else
      @ip_counter = (ips.size) - 1   # - c.c.c.c
      userinfo  = Userinfo.find(:first, :conditions => {:ip => "c.c.c.c"})
      @submission_counter = userinfo.overall_submissions
    end
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
    @speed = params[:speed][:speed]
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

    @job_description = params[:job_desc].gsub(/\s/,"__"); ### save that nobody enters sql syntax
    
    buildJobDir
    @raxml = Raxml.new({ :alifile =>params[:raxml][:alifile] , :query => @query, :outfile => @outfile, :speed => @speed, :substmodel => @substmodel, :heuristic => @heuristic, :treefile => params[:treefile][:file], :email => @email, :h_value => @h_value, :errorfile => "", :use_heuristic => @use_heuristic, :use_bootstrap => @use_bootstrap, :b_random_seed => @b_random_seed, :b_runs => @b_runs , :parfile => @parfile, :use_queryfile => @use_queryfile, :queryfile => @queryfile, :use_clustering => @use_clustering, :jobid => @jobid, :user_ip => @ip, :job_description => @job_description, :status => "running"})
    
    
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
      render :action => 'submit'
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
    @ip_counter = 0;
    @submission_counter = 0;
    getInfo
    @raxml = Raxml.find(:first, :conditions => ["jobid = #{params[:id]}"])
    @id = params[:id]
    if !(jobIsFinished?(@raxml.jobid))
      render :action => "wait"
    else
      @raxml.update_attribute(:status,"done")
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
    @cites = []
    jobid = params[:id]
    collectCites(jobid)
    @ip_counter = 0;
    @submission_counter = 0;
    getInfo
    rax =  Raxml.find(:first, :conditions => ["jobid = #{jobid}"])
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

  def collectCites(jobid)
    @cites << "<b>EPA:</b> <li> S.A. Berger, A. Stamatakis, Evolutionary Placement of Short Sequence Reads. arXiv:0911.2852v1 [q-bio.GN](2009)</li>"
    @cites << "<b>Archaeopteryx Treeviewer:</b> <li>Han, Mira V.; Zmasek, Christian M. (2009). phyloXML: XML for evolutionary biology and comparative genomics. BMC Bioinformatics (United Kingdom: BioMed Central) 10: 356. doi:10.1186/1471-2105-10-356. http://www.biomedcentral.com/1471-2105/10/356.</li>"
    @cites << "<li>Zmasek, Christian M.; Eddy, Sean R. (2001). ATV: display and manipulation of annotated phylogenetic trees. Bioinformatics (United Kingdom: Oxford Journals) 17 (4): 383–384. http://bioinformatics.oxfordjournals.org/cgi/reprint/17/4/383.</li>"
    rax =  Raxml.find(:first, :conditions => ["jobid = #{jobid}"])
    if rax.use_clustering.eql?("T")
      @cites << "<b>Hmmer:</b> <li>S. R. Eddy., A New Generation of Homology Search Tools Based on Probabilistic Inference. Genome Inform., 23:205-211, 2009.</li>"
      @cites << "<b>uclust:</b> <li>http://www.drive5.com/uclust</li>"
    end

  end    
  def download 
    file = params[:file]
    send_file file
  end

  def index
    getInfo
  end
  
  def look
    getInfo
    @error_id = ""
    @error_email = ""
    if !(params[:id].nil?)
      @error_id = "The job  with the id \'#{params[:id]}\' does not exists or is not finished yet"
    elsif !(params[:email].nil?) && !(params[:email].eql?("\'\'"))
      @error_email = "No jobs for \'#{params[:email]}\' available!"
    end
  end

  def findJob
    jobid = params[:rax_job]
    if Raxml.exists?(:jobid => jobid)
      if jobIsFinished?(jobid)
        redirect_to :action => "results" , :id => jobid
      else
        redirect_to :action => "look" ,:id => jobid
      end
    else
      redirect_to :action => "look" ,:id => jobid
    end
  end

  def listOldJobs
    jobs_email = params[:rax_email]
    if Raxml.exists?(:email => jobs_email) && (!jobs_email.eql?(""))
     
      redirect_to :action => "allJobs" , :email =>  "\'#{jobs_email}\'"
    else
      redirect_to :action => "look" ,:email => "\'#{jobs_email}\'"
    end
  end

  def allJobs
    getInfo
    jobs_email = params[:email]
    rax =  Raxml.find(:all, :conditions => ["email = #{jobs_email}"])
    @jobids=[]
    @jobdescs=[]
    rax.each do |r| 
      if jobIsFinished?(r.jobid)
        if r.job_description.eql?("") 
          @jobids << r.jobid
          @jobdescs << "";
        else
          @jobids << r.jobid
          @jobdescs << "<br/>"+r.job_description.gsub(/__/," ")
        end
      end
    end
  end
  
  def contact
    getInfo
    @error = ""
    if !(params[:id].nil?)
      @error = "An error occurres, please try again!"
    end
  end

  def sendMessage
    name = params[:con_name]
    name = name.gsub(/\s/,"__")
    email = params[:con_email]
    subject = params[:con_subject]
    subject = subject.gsub(/\s/,"__")
    subject = subject.gsub(/\"/,"\\\"")
    subject = subject.gsub(/\'/,"\\\\\'")
    message = params[:con_message]
    message = message.gsub(/\n/,"#n#")
    message = message.gsub(/\s/,"__")
    message = message.gsub(/\"/,"\\\"")
    message = message.gsub(/\'/,"\\\\\'")
    if Raxml.sendMessage(name,email,subject,message)
      redirect_to :action => "confirmation"
    else
      redirect_to :action => "contact", :id=>1
    end
  end

  def confirmation
    getInfo
  end

  def about
    getInfo
  end

 
end
