
require 'pty'
#RAILS_ROOT = "/home/denis/work/RAxMLWS"

class RaxmlTreefileParser
  
attr_reader :format, :valid_format, :error ,:data
  def initialize(stream)
   # @filename = stream.original_filename
    @filename ="file"
    @data = stream.readlines 
    @format = "unk"
    @valid_format = false
    @error = ""
    @message = "Invalid Treefile format!\n ParserError :: #{@filename}"
    check_format   
  end

  private
  def check_format
    @valid_format = true
    @format = "tree"
    checkTreefileFormatWithJava
  end

  private
  def checkTreefileFormatWithJava
    
    random_number = (1+rand(10000))* (1+(10000%(1+rand(10000))))*(1+rand(10000)) #build random number for @filename to avoid collision
    file = "#{RAILS_ROOT}/tmp/files/#{random_number}_#{@filename}" 
    f = File.open(file,'wb')
    @data.each {|d| f.write(d)}
    f.close
    cmd = "java -jar #{RAILS_ROOT}/bioprogs/java/treecheck.jar #{file}"
    # let RAxML check if phylip format is correct
    PTY.spawn(cmd) do |stdin, stdout, pid| 
      
      stdin.each do  |line| 
        if !(line =~ /good/)
          @error = " "
          @format = "unk"
          @valid_format = false
        end
      end
    end rescue Errno::EIO
    if !@error.eql?("")
      @error = "#{@message}\n#{@error}"
    end
    system "rm #{file}"
  end
end
