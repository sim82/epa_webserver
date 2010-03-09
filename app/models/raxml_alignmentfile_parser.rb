
require 'pty'

class RaxmlAlignmentfileParser
  
attr_reader :format, :valid_format, :error ,:data
  def initialize(stream)
    @filename = stream.original_filename
    @data = stream.readlines 
    @format = "unk"
    @valid_format = false
    @error = ""
    check_format   
  end

  private
  def check_format
    i = 0
    while i < @data.size

      if @data[i] =~ /^>/ && @data[i+1]=~ /^[A-Z\-]+/ #fasta
        @format = "fas"
        @valid_format = true
        j = i+2
        while j < @data.size
          if !(@data[j] =~ /^[A-Z\-]+/ || (@data[j] =~  /^>/ && @data[j+1]=~ /^[A-Z\-]+/))
            @format = "unk"
            @valid_format = false
            @error = "Unknown file format!(File seems to be fasta formatted)\n ParserError :: #{@filename} line: #{j}\n *#{@data[j]}*"
            return
          end
        j = j+1
        end

      elsif @data[i] =~ /^\d+\s\d+/  #phylip
        @format = "phl"
        @valid_format = true
        random_number = (1+rand(10000))* (1+(10000%(1+rand(10000))))*(1+rand(10000)) #build random number for @filename to avoid collision
        file = "#{RAILS_ROOT}/tmp/files/#{random_number}_#{@filename}" ### change before production phase
        f = File.open(file,'wb')
        @data.each {|d| f.write(d)}
        f.close
        cmd = "#{RAILS_ROOT}/bioprogs/raxml/raxmlHPC  -s #{file} -fc -m GTRGAMMA -n test"
        # let RAxML check if phylip format is correct
        PTY.spawn(cmd) do |stdin, stdout, pid| 
	  
          stdin.each do  |line| 
            if !(line =~ /^Alignment\sformat\scan\sbe\sread\sby\sRAxML/)
              @error = @error+line
              @format = "unk"
              @valid_format = false
            end
          end
        end rescue Errno::EIO
        if !@error.eql?("")
           @error = "Unknown file format!(File seems to be phylip formatted)\n ParserError :: #{@filename}\n#{@error}"
        end
        system "rm #{file}"
        break

      elsif @data[i] =~ /^\s+$/ # blank lines are ignored
        #do nothing

      else
        @error = "Unknown format, file not parsable!"
        break
      end
      i = i+1
    end
  end
end

