function readsettings()
    str = open(f -> read(f, String), EmailPass_File_path)
    settings = JSON.parse(str)
    return settings
end


function sendmail(filename::String, scf_po::Bool, dft_mixing::Mixing)

    println("")
    println("Sending results")

    if scf_po
        subject = "(SCF) $(filename) calculation convergence"
    else
        subject = "(SCF) $(filename) calculation not convergence"
    end

    message = "SCF results"
        
    make_attachments(dft_mixing)
    sendmail(subject, message, ["temp_attachment.txt"])
end


function sendmail(filename::String, GeoOpt_po::Bool, geo_optim::Geo_Optim)

    println("")
    println("Sending results")

    if GeoOpt_po
        subject = "(Opt) $(filename) calculation convergence"
    else
        subject = "(Opt) $(filename) calculation not convergence"
    end

    message = "Geometric Optimization results"
        
    make_attachments(geo_optim)
    sendmail(subject, message, ["temp_attachment.txt"])
end


# For attachments
# sendmail("test send .std file", "test send results", ["Cdia.txt"])
function sendmail(subject::String, message::String, attachments::Vector{String})
    settings = readsettings()
    opt = SendOptions(
        isSSL=true,
        username=settings["username"],
        passwd=settings["passwd"])
    # Provide the message body as RFC5322 within an IO
    dtstr = Dates.format(Dates.now(), "e, dd u YYYY HH:MM:SS +0900")
    body = IOBuffer(
        "Date: " * dtstr * "\r\n" *
        "From: " * settings["from"] * " <" * settings["username"] * ">\r\n" *
        "To: " * settings["username"] * "\r\n" *
        "Subject: " * subject * "\r\n" *
        "\r\n" *
        message * 
        "\r\n")
    url = settings["url"]
    rcpt = ["<" * settings["username"] * ">"]
    to = ["<" * settings["username"] * ">"]
    from = "<" * settings["username"] * ">"

    mime_msg = get_mime_msg(message)
    body = get_body(to, from, subject, mime_msg; attachments)
    resp = send(url, rcpt, from, body, opt)
end

