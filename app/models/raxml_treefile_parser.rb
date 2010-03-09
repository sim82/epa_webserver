class RaxmlTreefileParser
  
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
        j = i+1
        while j < @data.size
          if !(@data[j] =~ /^\w+\d+\s+[A-Z\-]+/)
            @format = "unk"
            @valid_format = false
            @error = "Unknown file format!(File seems to be phylip formatted)\n ParserError :: #{@filename} line: #{j}\n *#{@data[j]}*"
            return
          end
        j = j+1
        end
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
