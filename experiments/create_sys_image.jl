using PackageCompiler
ENV["PROFILE"] = "1"
PackageCompiler.create_sysimage(["OceanScalingTests"]; sysimage_path="OST-sysimage.so", precompile_statements_file="run.jl")

