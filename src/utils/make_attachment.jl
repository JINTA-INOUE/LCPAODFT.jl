function make_attachments(dft_mixing::Mixing)

    is_convergence = dft_mixing.is_convergence
    Conv_iter = dft_mixing.Conv_iter
    NormRD = dft_mixing.NormRD
    HisEele = dft_mixing.HisEele

    data = open(pwd()*"/"*PROGRAM_FILE, "r")
    julia_inputfile = readlines(data)
    close(data)

    attachments = open("temp_attachment.txt", "w")
    println(attachments, julia_inputfile)
    println(attachments, "\n\n\n")
    println(attachments, "============== SCF results ============== ")
    println(attachments, "\t\t is_convergence = $is_convergence")
    println(attachments, "\t\t Conv_iter = $Conv_iter")
    println(attachments, "\n")
    for iter = 1:Conv_iter
        @printf(attachments, "iter = %3d  NormRD = %5.12f  Eele = %5.12f", iter, NormRD[iter], HisEele[iter])
    end
    println(attachments, "========================================= ")

    close(attachments)
end

#=
function make_attachments(geo_optim::Geo_Optim)

    Natom = geo_optim.Natom
    Geo_Opt_convergence = geo_optim.Geo_Opt_convergence
    Geo_Opt_Conv_iter = geo_optim.Geo_Opt_Conv_iter
    GxyzHisIn = geo_optim.GxyzHisIn
    GxyzHisR = geo_optim.GxyzHisR
    Geo_Opt_Max_Force = geo_optim.Geo_Opt_Max_Force

    data = open(pwd()*"/"*PROGRAM_FILE, "r")
    julia_inputfile = readlines(data)
    close(data)

    attachments = open("temp_attachment.txt", "w")
    println(attachments, julia_inputfile)
    println(attachments, "\n\n\n")
    println(attachments, "============== Geometric Optimization results ============== ")
    println(attachments, "\t\t Geo_Opt_convergence = $Geo_Opt_convergence")
    println(attachments, "\t\t Geo_Opt_Conv_iter = $Geo_Opt_Conv_iter")
    println(attachments, "")
    println(attachments, "<Max_Force>")
    for iter = 1:Geo_Opt_Conv_iter
        @printf(attachments, "iter = %3d  Max_Force = %5.12f\n", iter, Geo_Opt_Max_Force[iter])
    end

    println(attachments, "")
    println(attachments, "<Gxyz>")
    for iter = 1:Geo_Opt_Conv_iter, atom = 1:Natom
        @printf(attachments, "iter = %3d  atom = %d   %5.12f  %5.12f  %5.12f\n", iter, GxyzHisIn[iter][atom][1], GxyzHisIn[iter][atom][2], GxyzHisIn[iter][atom][3])
    end

    println(attachments, "")
    println(attachments, "<Force>")
    for iter = 1:Geo_Opt_Conv_iter, atom = 1:Natom
        @printf(attachments, "iter = %3d  atom = %d   %5.12f  %5.12f  %5.12f\n", iter, GxyzHisR[iter][atom][1], GxyzHisR[iter][atom][2], GxyzHisR[iter][atom][3])
    end
    println(attachments, "============================================================= ")

    close(attachments)
end
=#