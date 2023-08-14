;; To get a development container (inside emacs shell will work)
;;
;;   guix shell -C -D -f guix.scm -- bash --init-file <(echo "ln -s /bin/sh /bin/bash")
;;
;; and build
;;
;;   find -name CMakeCache.txt|xargs rm -v
;;   cd build
;;   cmake -DCMAKE_BUILD_TYPE=Debug ..
;;   cmake --build . --verbose -- -j 14 && ctest . --verbose
;;
;; For the tests you may need /usr/bin/env. In a container create it with
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;

(use-modules
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim)
  ((guix licenses) #:prefix license:)
  (guix packages)
  (guix utils)
  (guix download)
  (guix git-download)
  (guix build-system gnu)
  (guix gexp)
  (gnu packages base)
  (gnu packages bioinformatics) ; for samtools in sambamba
  (gnu packages compression)
  (gnu packages curl)
  (gnu packages dlang)
  (gnu packages gcc)
  (gnu packages perl)
  (gnu packages python)
  (gnu packages ruby)
  (gnu packages tls)
  (gnu packages version-control)
  )

#!
 (use-modules
  (guix build-system cmake)
  (guix utils)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages build-tools)
  (gnu packages commencement) ; gcc-toolchain
  (gnu packages curl)
  (gnu packages datastructures)
  (gnu packages gdb)
  (gnu packages gcc)
  (gnu packages jemalloc)
  (gnu packages libffi)
  (gnu packages mpi)
  (gnu packages python)
  (gnu packages python-xyz)
  (gnu packages pkg-config)
  (gnu packages tls)
)

!#

(define %source-dir (dirname (current-filename)))

(define %git-commit
  (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public sambamba-git
  (package
    (name "sambamba-git")
    (version (git-version "1.0.0" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
      (build-system gnu-build-system)
      (outputs '("out"     ; disable all checks for speed
                 "debug"))
      (inputs
       `(("samtools" ,samtools) ; for pileup
         ("bcftools" ,bcftools) ; for pileup
         ("lz4" ,lz4)
         ("zlib" ,zlib)
       ))
      (native-inputs
       `(("ldc" ,ldc)
         ("coreutils" ,coreutils) ; for env
         ("perl" ,perl) ; Needed for building htslib
         ("ruby" ,ruby) ; Needed for building htslib
         ("python" ,python-2) ; Needed for building htslib and sambamba
         ("gcc" ,gcc)
         ("which" ,which)))
      (home-page "https://github.com/lomereiter/sambamba")
      (synopsis "Fast tool for working with SAM, BAM and CRAM files written in D.")
      (description
       "Sambamba is a high performance modern robust and fast
tool (and library), written in the D programming language, for working
with SAM, BAM and CRAM files.  Current parallelised functionality is
an important subset of samtools functionality, including view, index,
sort, markdup, and depth.")
      (license license:gpl2+)))

sambamba-git
