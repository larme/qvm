;;; package.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(defpackage #:dqvm2
  (:use #:common-lisp
        #:cl-mpi
        #:cl-mpi-extensions
        #:static-vectors)
  (:local-nicknames (#:a #:alexandria))
  (:import-from #:qvm
                #:boolean-bit
                #:defun-inlinable
                #:number-of-qubits
                #:transition)
  (:export #:distributed-qvm
           #:make-distributed-qvm
           #:print-qubit-permutation
           #:apply-qubit-permutation
           #:apply-inverse-qubit-permutation
           #:rank
           #:number-of-processes
           #:number-of-qubits
           #:permutation
           #:addresses
           #:number-of-amplitudes
           #:reset-wavefunction
           #:save-wavefunction
           #:find-next-addresses
           #:update-permutation
           #:qubit-permutation
           #:number-of-addresses
           #:do-addresses
           #:address-memberp
           #:offset-of
           #:get-rank
           #:get-address-by-offset
           ))

(defpackage #:dqvm2-user
  (:use #:common-lisp
        #:cl-mpi
        #:cl-mpi-extensions
        #:dqvm2))
