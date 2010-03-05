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
    heuristics = ["none","MP","MC"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    @heuristics_values =""
    heuristics_values = ["0.5","0.333","0.25","0.2","0.167","0.142","0.125","0.111","0.1"]
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
    heuristics = ["none","MP","ML"]
    heuristics.each {|h| @heuristics = @heuristics+"<option>#{h}</option>"}
    @heuristics_values =""
    heuristics_values = ["0.5","0.333","0.25","0.2","0.167","0.142","0.125","0.111","0.1"]
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
      @substmodel = "#{@substmodel}_#{@matrix}_#{@sm_float}"
    end
    
    @heuristic = params[:heuristic]
    @heu_float = "" 
    if !@heuristic.eql?("none")
      @heu_float = params[:heu_float]
      @heuristic = @heuristic+"_"+@heu_float
    end
#    @heuristic = @heuristic+"_"+@heu_float

    @wait = params[:wait]
    @email = params[:rax_email]
    @outfile = ""
    @alifile = ""
    @pid = 0

    @raxml = Raxml.new({ :alifile =>params[:raxml][:alifile] , :query => @query, :outfile => @outfile, :speed => @speed, :substmodel => @substmodel, :heuristic => @heuristic, :treefile => params[:treefile][:file], :email => @email, :pid => @pid, :wait => @wait})
    
    
    if @raxml.save
      buildJobDir(@raxml)
      @alifile = saveInfile(params[:raxml][:alifile])
      @raxml.update_attribute(:alifile,@alifile)
      @raxml.update_attribute(:outfile,@directory+"results.txt")

      @treefile = saveInfile(params[:treefile][:file])
      @raxml.update_attribute(:treefile,@treefile)
      
      link = url_for :controller => 'raxml', :action => 'results', :id => @raxml.id
      @raxml.execude(link)
      sleep 2
      redirect_to :action => 'wait', :id => @raxml.id 
    else
      puts "#############"
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

  def saveInfile(stream)
    file = @directory+stream.original_filename
    File.open(file, "wb") { |f| f.write(stream.read) }
    return file
  end
   
  def wait
    rax = Raxml.find_by_id(params[:id])
    if pidAlive?(rax.pid)
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

  def results
    rax = Raxml.find_by_id(params[:id])    
#    if rax.emailValid?
#      sendEmail(rax.email,rax.id)
#    end
    @results =  RaxmlResultsParser.new(rax.outfile).data
  end

  def sendEmail(recipient,id) ### not in use
    subject = "Your RAxML job has been finished"
    link = url_for :controller => 'raxml', :action => 'results', :id => id
    Emailer.deliver_email(recipient, subject, link)
    return if request.xhr?
  end
  
  def saveEmail ### not in use
    rax = Raxml.find_by_id(params[:id])
    rax.update_attribute(:email, params[:email][:e])
    redirect_to :action => 'wait' , :id => params[:id]
  end
end
