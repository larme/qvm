(defpackage #:dqvm2-system
  (:use #:common-lisp #:asdf))

(in-package #:dqvm2-system)

(defsystem #:dqvm2
  :description "Rigetti Distributed Quantum Virtual Machine"
  :author "Juan M. Bello Rivas <jbellorivas@rigetti.com>, Robert Smith <robert@rigetti.com>, Lauren Capelluto <lauren@rigetti.com>"
  :licence "Apache License 2.0 (See LICENSE.txt)"
  :depends-on (
               #:alexandria
               #:cl-mpi
               #:cl-mpi-extensions
               #:cl-quil
               #:cl-syslog
               #:command-line-arguments
               #:magicl
               #:qvm
               ;; #:regex-trace
               #:static-vectors
               #:trivial-garbage
               )
  :defsystem-depends-on (#:cl-mpi-asdf-integration)
  :class :mpi-program
  :build-operation :static-program-op
  :build-pathname "dqvm2"
  :entry-point "dqvm2:main"
  :serial t
  :pathname "src/"
  :components ((:file "package")
               (:file "utils")
               (:file "mpi")
               (:file "logging")
               (:file "addresses")
               (:file "distributed-qvm")
               (:file "linear-algebra")
               (:file "permutations")
               (:file "sendrecv")
               (:file "measurement")
               (:file "apply-distributed-gate")
               (:file "transition")
               (:file "entry-point"))
  :in-order-to ((test-op (test-op "dqvm2/tests"))))

(defsystem #:dqvm2/tests
  :description "Test suite for dqvm2."
  :depends-on (#:dqvm2 #:fiasco)
  :components
  ((:file "test-suite")))

(defmethod perform ((op test-op) (system (eql (find-system "dqvm2"))))
  (let ((run-tests (find-symbol "ALL-TESTS" "FIASCO")))
    (funcall run-tests)))

(defmethod operation-done-p ((op test-op) (system (eql (find-system "dqvm2"))))
  nil)

(delete-package *package*)
