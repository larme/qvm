;;; distributed-qvm.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package #:dqvm2)

(defclass distributed-qvm (qvm:classical-memory-mixin)
  ((addresses
    :accessor addresses
    :initarg :addresses
    :initform (error-missing-initform :addresses)
    :type addresses
    :documentation "Address table.")
   (amplitudes
    :accessor amplitudes
    :initform nil
    :type static-vector-pointer        ; qvm::quantum-state
    :documentation "The components of the (permuted) wave function.")
   (scratch
    :accessor scratch
    :initform nil
    :type static-vector-pointer
    :documentation "Temporary memory to store amplitudes."))

  (:documentation "A distributed implementation of the Quantum Abstract Machine. A DQVM object keeps track of data owned by a single MPI rank."))

(defun make-distributed-qvm (&rest rest)
  "Convenience function for instantiating DISTRIBUTED-QVM objects. The parameters are the same as the intargs of ADDRESSES."
  (let* ((addresses (apply #'make-instance 'addresses rest)))
    (make-instance 'distributed-qvm :addresses addresses)))

(defmethod initialize-instance :after ((qvm distributed-qvm) &rest initargs)
  (declare (ignore initargs))

  (reset-wavefunction qvm)

  (let ((scratch (scratch qvm))
        (amplitudes (amplitudes qvm)))
    (trivial-garbage:finalize qvm (lambda ()
                                    (when scratch
                                      (free-static-vector scratch))
                                    (when amplitudes
                                      (free-static-vector amplitudes))))))

(defmethod qvm:number-of-qubits ((qvm distributed-qvm))
  (number-of-qubits (addresses qvm)))

(defmethod rank ((qvm distributed-qvm))
  (rank (addresses qvm)))

(defmethod number-of-processes ((qvm distributed-qvm))
  (number-of-processes (addresses qvm)))

(defun make-complex-static-vector (n)
  "Create empty complex static vector of size N."
  (make-static-vector n :element-type 'qvm:cflonum
                        :initial-element (qvm:cflonum 0)))

(defmethod reset-wavefunction ((qvm distributed-qvm))
  (with-slots (amplitudes scratch) qvm

    (let* ((rank (rank qvm))
           (number-of-processes (number-of-processes qvm))
           (number-of-qubits (number-of-qubits qvm))
           (maximum-arity (get-maximum-arity (qvm::program qvm)))
           (block-size (if (plusp maximum-arity)
                           (* 2 maximum-arity)
                           *default-block-size*))
           (addresses (make-instance 'addresses :rank rank
                                                :number-of-processes number-of-processes
                                                :number-of-qubits number-of-qubits
                                                :block-size block-size))
           (number-of-addresses (number-of-addresses addresses)))

      ;; Initialize address table.
      (setf (addresses qvm) addresses)

      ;; Initialize wave function proper.
      (when amplitudes
        (free-static-vector amplitudes))
      (setf amplitudes (make-complex-static-vector number-of-addresses))
      (when (address-memberp addresses 0)
        (setf (aref (amplitudes qvm) (offset-of addresses 0))
              (qvm:cflonum 1)))

      ;; Initialize temporary memory.
      (when scratch
        (free-static-vector scratch))
      (setf scratch (make-complex-static-vector number-of-addresses))))

  qvm)

(defmethod reset-wavefunction-debug ((qvm distributed-qvm))
  (with-slots (amplitudes addresses)
      qvm
    (do-addresses (address addresses)
      (setf (aref amplitudes (offset-of addresses address))
            (qvm:cflonum address)))))

(defun load-instructions (qvm instructions)
  "Load the code vector INSTRUCTIONS into QVM."
  ;; XXX ensure the qubits needed match the current qvm qubit count.
  (setf (qvm::program qvm) instructions)
  (reset-wavefunction qvm)
  qvm)

(defmethod print-object ((qvm distributed-qvm) stream)
  (let ((*print-readably* nil)
        (*print-pretty* nil))
    (print-unreadable-object (qvm stream :type t :identity t)
      (with-slots (addresses amplitudes scratch)
          qvm
        (format stream "~{~A ~A~^ ~}"
                (list (prin1-to-string :addresses) addresses
                      (prin1-to-string :amplitudes) amplitudes
                      (prin1-to-string :scratch) scratch))))))

(defun qubit-permutation (instruction)
  "Return an association list representing the qubits we must transpose to execute INSTRUCTION."
  (check-type instruction quil:instruction)
  (loop :for i :from 0 :for qubit :in (quil::arguments instruction)
        :when (/= i (quil:qubit-index qubit))
        :collect (list i (quil:qubit-index qubit))))

(defmethod save-wavefunction ((qvm distributed-qvm) filename)
  "Save wavefunction in QVM to FILENAME in parallel.

The file format is a consecutive set of bytes containing real and imaginary parts of the (ordered) amplitudes."
  (format-log :debug "Saving wave function to ~S" filename)

  (uiop:delete-file-if-exists filename)

  (with-slots (addresses amplitudes)
      qvm

    (let ((amode (logior mpi::+mpi-mode-create+ mpi::+mpi-mode-wronly+)))

      (mpi::with-open-mpi-file (fh filename amode)
        (loop :for offset :from 0 :by qvm::+octets-per-cflonum+
              :for ptr-amplitude := (static-vector-pointer amplitudes :offset offset)
              :for index :from 0 :below (number-of-addresses addresses)
              :for address := (get-address-by-offset addresses index)
              :for file-offset := (* address qvm::+octets-per-cflonum+)
              :do (mpi::%mpi-file-write-at fh file-offset ptr-amplitude 1 +mpi-cflonum+
                                           (cffi:null-pointer)))))))

;;; XXX Add functions to handle the case when a program is reloaded and block-size changes.
