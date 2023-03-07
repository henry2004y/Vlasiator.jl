# Vlasiator automatic installation

using ArgParse

function parse_commandline()
   s = ArgParseSettings()

   @add_arg_table s begin
      "--install", "-i"
         help = "install Vlasiator"
         action = :store_true
      "--installdep", "--idep"
         help = "install Vlasiator dependencies"
         action = :store_true
      "--lib"
         help = "install Vlasiator dependency name"
         arg_type = Symbol
         default = :all
      "--dir"
         help = "installation location"
         arg_type = String
         default = "lib"
   end

   return parse_args(s)
end

function main()
   parsed_args = parse_commandline()

   if parsed_args["install"]
      parsed_args["installdep"] = true
   end

   # Install dependencies
   if parsed_args["installdep"]
      !isdir(parsed_args["dir"]) && mkdir(parsed_args["dir"])
      cd(parsed_args["dir"])
      origindir = pwd()
      tmpdir = mktempdir(origindir)
      cd(tmpdir)

      if parsed_args["lib"] == :all || parsed_args["lib"] == :dccrg
         libpath = joinpath(origindir, "dccrg")
         if !isdir(libpath)
            run(`git clone https://github.com/fmihpc/dccrg.git`)
            cd("dccrg")
            run(`git checkout vlasiator-version`)
            cd("../")
            mv("dccrg", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :vlsv
         libpath = joinpath(origindir, "vlsv")
         if !isdir(libpath)
            run(`git clone https://github.com/fmihpc/vlsv.git`)
            cd("vlsv")
            run(`make`)
            cd("../")
            mv("vlsv", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :phiprof
         libpath = joinpath(origindir, "phiprof")
         if !isdir(libpath)
            run(`git clone https://github.com/fmihpc/phiprof.git`)
            cd("phiprof/src")
            run(`make`)
            cd("../../")
            mv("phiprof", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :vlsv
         libpath = joinpath(origindir, "vlsv")
         if !isdir(libpath)
            run(`git clone https://github.com/fmihpc/vlsv.git`)
            cd("vlsv")
            run(`make`)
            cd("../")
            mv("vlsv", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :vectorclass
         joinpath(origindir, "vectorclass")
         if !isdir(libpath)
            run(`git clone https://github.com/vectorclass/version1.git`)
            run(`git clone https://github.com/vectorclass/add-on.git`)
            cp("add-on/vector3d/vector3d.h", "version1/vector3d.h")
            mv("version1", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :jemalloc
         libpath = joinpath(origindir, "jemalloc")
         if !isdir(libpath)
            run(`wget https://github.com/jemalloc/jemalloc/releases/download/4.0.4/jemalloc-4.0.4.tar.bz2`)
            run(`tar -xf jemalloc-4.0.4.tar.bz2`)
            cd("jemalloc-4.0.4")
            run(`./configure --prefix=$libpath --with-jemalloc-prefix=je_`)
            run(`make`)
            run(`make install`)
            cd("../")
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :eigen
         libpath = joinpath(origindir, "Eigen")
         if !isdir(libpath)
            run(`wget https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.bz2`)
            run(`tar -xf eigen-3.2.8.tar.bz2`)
            mv("eigen-3.2.8/Eigen", libpath)
         end
      end

      if parsed_args["lib"] == :all || parsed_args["lib"] == :zoltan
         libpath = joinpath(origindir, "zoltan")
         if !isdir(libpath)
            run(`git clone https://github.com/sandialabs/Zoltan.git`)
            mkdir(libpath)
            cd(libpath)
            zoltanpath = joinpath(tmpdir, "Zoltan")
            run(`$zoltanpath/configure --prefix=$libpath --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong`)
            run(`make`)
            run(`make install`)
         end
      end
   end

   if parsed_args["install"]
      run(`git clone https://github.com/fmihpc/vlasiator.git`)
   end
end

main()