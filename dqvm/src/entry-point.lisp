;;; entry-point.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package #:dqvm2)

(defun %main (argv)                       ; XXX change name, follow qvm-app.
  (when (zerop +rank+)
    (format-log :debug "Command line arguments: \"~{~A~^ ~}\"" argv)

    (unless (rest argv)
      (format uiop/stream:*stderr* "Usage: ~A QUIL-FILE~%" (first argv))
      (uiop:quit -1)))

  (let (qvm)

    (cond
      ((zerop +rank+)
       (let ((filename (second argv)))
         (format-log :info "Reading Quil file ~S" filename)
         (let* ((program (quil:read-quil-file filename))
                (number-of-qubits (quil:qubits-needed program)))
           (format-log :info "Read program ~S using ~D total qubits. Code: ~{~A~^, ~}"
                       filename number-of-qubits
                       (map 'list #'instruction->string
                            (quil:parsed-program-executable-code program)))
           (bcast :value number-of-qubits)
           (setf qvm (make-distributed-qvm :number-of-qubits number-of-qubits))
           (load-instructions qvm (bcast-instructions program)))))
      (t
       (format-log :info "Waiting for instructions")
       (let ((number-of-qubits (bcast))
             (instructions (bcast-instructions)))
         (setf qvm (make-distributed-qvm :number-of-qubits number-of-qubits))
         (load-instructions qvm instructions))))

    ;; (reset-wavefunction-debug qvm)

    (let ((seed (+ (nth-value 1 (sb-ext:get-time-of-day)) ; XXX SBCLism.
                   +rank+)))
      (qvm:with-random-state ((qvm:seeded-random-state seed))
        (qvm:run qvm)))

    (save-wavefunction qvm "wavefunction.dat")
    (print-reference-result qvm (second argv))

    (format-log :info "Finished program execution.")))

(defun print-reference-result (qvm filename)
  (when (zerop +rank+)
    (let* ((qvm:*transition-verbose* nil)
           (matrix (qvm:program-matrix (quil:read-quil-file filename)))
           (m (magicl:matrix-cols matrix))
           (components (loop :for i :from 0 :below m
                             ;; :collect (qvm:cflonum i)
                             :collect (if (zerop i)
                                          (qvm:cflonum 1)
                                          (qvm:cflonum 0))))
           (initial-wavefunction (magicl:make-complex-matrix m 1 components))
           (wavefunction (magicl:multiply-complex-matrices matrix initial-wavefunction)))

      (format-log :info "Expected wavefunction:~%~{~A~%~}"
                  (coerce (magicl::matrix-data wavefunction) 'list)))))

(defun entry-point (argv)
  (let ((*print-pretty* nil)
        (*print-case* :downcase)
        (qvm:*transition-verbose* t) ; XXX turn this on/off via command line args.
        )

    (uiop/stream:setup-stderr)

    ;; (trace compute-matrix-vector-products)
    ;; (regex-trace:regex-trace "^%MPI-.?(SEND|RECV|WAIT)" :print (mpi-comm-rank))

    (unless (mpi-initialized)
      (mpi-init :thread-support :mpi-thread-multiple))

    (mpi::%mpi-comm-set-errhandler +mpi-comm-world+ +mpi-errors-are-fatal+)

    (setf +rank+ (mpi-comm-rank))

    (setup-logger "Welcome to the Rigetti Distributed Quantum Virtual Machine")
    (let ((status 0))
      (unwind-protect
           (setf status (%main argv))
        (mpi-finalize)
        (uiop:quit status)))))
