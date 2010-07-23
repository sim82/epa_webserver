
require 'pty'
require 'fasta_to_phylip'

class RaxmlAlignmentfileParser
  
attr_reader :format, :valid_format, :error ,:data , :ali_length
  def initialize(stream)
    @filename = stream.original_filename
    @data = stream.readlines 
    @format = "unk"
    @valid_format = false
    @error = ""
    @message = "Invalid alignmentfile format!(Only Fasta and Phylip alignment formats are allowed)  ParserError :: #{@filename}<br></br>"
    check_format   
    @ali_length = 0
    if @valid_format
      @data.each do |line|
        if line=~ /\d+\s+(\d+)/
          @ali_length = $1.to_i
        end
      end
    end
  end

  private
  def check_format
    i = 0
    while i < @data.size

      if @data[i] =~ /^>/ && @data[i+1]=~ /^([A-Z\-]+)/ #fasta format?
        @format = "fas"
        @valid_format = true
        j = i+2
        while j < @data.size
          if !(@data[j] =~ /^[A-Z\-]+/ || (@data[j] =~  /^>/ && @data[j+1]=~ /^[A-Z\-]+/))
            @format = "unk"
            @valid_format = false
            @error = "#{message} line: #{j}\n *#{@data[j]}*"
            return
          end
        j = j+1
        end
        @data = FastaToPhylip.new(@data).phylip.join("\n")
        checkPhylipFormatWithRaxml
        break

      elsif @data[i] =~ /\s*\d+\s+\d+/  #phylip format?
        @format = "phl"
        @valid_format = true
        checkPhylipFormatWithRaxml
        break

      elsif @data[i] =~ /^\s+$/ # blank lines are ignored
        #do nothing

      else
        @error = @message
        break
      end
      i = i+1
    end
  end

  private
  def checkPhylipFormatWithRaxml
    
    random_number = (1+rand(10000))* (1+(10000%(1+rand(10000))))*(1+rand(10000)) #build random number for @filename to avoid collision
    file = "#{RAILS_ROOT}/tmp/files/#{random_number}_#{@filename}" 
    f = File.open(file,'wb')
    @data.each {|d| f.write(d)}
    f.close
    cmd = "#{RAILS_ROOT}/bioprogs/raxml/raxmlHPC  -s #{file} -fc -m GTRGAMMA -n test"
    #let RAxML check if phylip format is correct
    PTY.spawn(cmd) do |stdin, stdout, pid| 
      
      stdin.each do  |line| 
        if !(line =~ /^Alignment\sformat\scan\sbe\sread\sby\sRAxML/) || !(line=~ /IMPORTANT\s+WARNING/)
          @error = @error+line
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

