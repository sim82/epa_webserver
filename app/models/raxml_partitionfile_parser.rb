class RaxmlPartitionfileParser

  attr_reader :data, :valid_format, :error

  def initialize(stream, ali_length)
    @filename = stream.original_filename
    @data = stream.readlines
    @valid_format = false
    @error  = ""
    @len = ali_length
    @occ = Array.new(@len,0)
    check_format
  end

  private
  def check_format
    n = 1
    reads = []
    @data.each do |line|
      if line =~ /^\s*$/
        
      elsif  line =~ /^[A-Z]+,\s+\w+\s+=(\s+\d+-\d+,)*(\s+\d+-\d+)$/ ||  line =~ /^[A-Z]+,\s+\w+\s+=(\s+\d+-\d+\\\d+,)*(\s+\d+-\d+\\\d+)$/
        digits = line.scan(/\d+\-\d+/)
        digits.each do |dig|
          if dig =~/(\d+)\-(\d+)/
            if $1.to_i > @len || $2.to_i > @len
              @error = "Read position to long(Alignment length is #{@len})! \n ParserError :: #{@filename} line: #{n} => #{line}"
              return
            end
          end
        end
        if line =~ /^[A-Z]+,\s+\w+\s+=(\s+\d+-\d+,)*(\s+\d+-\d+)$/
          a = line.scan(/\d+\-\d+/)
          if line =~/(\d+)\-(\d+)/
            reads << [$1.to_i,$2.to_i]
          end
        end
      else
        @error = "Invalid partitionfile format! \n ParserError :: #{@filename} line: #{n} => #{line}"
        return
      end
      n = n+1
    end
    reads.each do |x,y|
      if !(occupy_sequence(x,y))
        @error = "Reads are not allowed to overlap! \n ParserError :: #{@filename}"
        return
      end
    end
    @valid_format = true
  end

  def occupy_sequence(i,j)
    i = i-1
    while i < j
      if @occ[i] == 1
        return false
      else
        @occ[i] = 1
      end
      i = i+1
    end
    return true
  end
end