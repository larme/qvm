;;; measurement.lisp
;;;
;;; Author: Juan M. Bello-Rivas
;;;         Robert Smith

(in-package #:dqvm2)

(defmethod measure ((qvm distributed-qvm) q)
  (check-type q qvm::nat-tuple-element)
  (assert (< q (number-of-qubits qvm)) (qvm q)
          "Trying to measure qubit ~D on a QVM with only ~D qubit~:P."
          q
          (qvm:number-of-qubits qvm))
  (let* ((r (mt19937:random 1.0d0))
         (excited-probability (wavefunction-excited-state-probability (addresses qvm) (amplitudes qvm) q))
         (cbit (mpi-broadcast-anything 0 :object (boolean-bit (<= r excited-probability)))))

    ;; Force the non-deterministic measurement.
    (force-measurement cbit q qvm excited-probability)

    ;; Return the qvm.
    (values qvm cbit)))

(defun wavefunction-excited-state-probability (addresses amplitudes q)
  (let ((number-of-addresses (number-of-addresses addresses)))
    (cffi:with-foreign-objects ((value :double)
                                (probability :double))
      (setf (cffi:mem-ref value :double)
            (loop :for offset :from 0 :below number-of-addresses ; XXX parallelize
                  :for address := (get-address-by-offset addresses offset)
                  :when (logbitp q address)
                    :summing (qvm:probability (aref amplitudes offset))
                  double-float))
      (mpi::%mpi-allreduce value probability 1 mpi:+mpi-double+ mpi:+mpi-sum+ mpi:*standard-communicator*)
      (cffi:mem-ref probability :double))))

(defun force-measurement (measured-value qubit qvm excited-probability)
  (let* ((wavefunction (amplitudes qvm))
         (annihilated-state (- 1 measured-value))
         (inv-norm (if (zerop annihilated-state)
                       (/ (sqrt excited-probability))
                       (/ (sqrt (- (qvm:flonum 1) excited-probability)))))
         (addresses (addresses qvm))
         (number-of-addresses (number-of-addresses addresses)))
    (loop :for offset :from 0 :below number-of-addresses ; XXX parallelize
          :for address := (get-address-by-offset addresses offset)
          :do (setf (aref wavefunction offset)
                    (if (= annihilated-state (ldb (byte 1 qubit) address))
                        (qvm:cflonum 0)
                        (* inv-norm (aref wavefunction offset))))))
  qvm)
