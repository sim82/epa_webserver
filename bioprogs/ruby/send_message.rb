#!/usr/bin/ruby

RAILS_ROOT = File.expand_path(File.join(File.dirname(__FILE__), '../..'))
require 'net/smtp'
require "#{File.dirname(__FILE__)}/../../config/environment.rb"
SERVER_NAME = ENV['SERVER_NAME']

class SendMessage 

  def initialize(opts)
    @message = opts.join(" ")
    @email_address1 = "Denis.Krompass@campus.lmu.de"
    @email_address2 = "stamatak@in.tum.de"
    @email_address3 = "bergers@in.tum.de"
    send_email
    
  end

  def send_email
    Net::SMTP.start('localhost', 25) do |smtp|
      smtp.open_message_stream("#{ENV['SERVER_NAME']}", @email_address1) do |f|
        
        f.puts "From: #{ENV['SERVER_NAME']}"
        
        f.puts "To: #{@email_address}"

        f.puts 'Subject: RAxMLWS message'
        
        f.puts "#{@message}"
      end
    end
  end
end

SendMessage.new(ARGV)
