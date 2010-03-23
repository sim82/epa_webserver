class RaxmlResultsParser

  attr_reader :data ,:names, :files

  def initialize(file_ending)
    @job_id = file_ending
    @data = []
    @names = []
    @files = []
    getFiles
  end

  def getFiles
    Dir.glob("#{RAILS_ROOT}/public/jobs/#{@job_id}/RAxML_*"){|file| 
      @files << file
      if file =~ /.+\/(RAxML_.+)\.#{@job_id}(\.*\d*)$/
        @names << $1+$2
      end
    }
  end

end
