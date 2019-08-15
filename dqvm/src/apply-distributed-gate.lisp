;;; apply-distributed-gate.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package :dqvm2)

(declaim (optimize (speed 3) (safety 0) (debug 0)))

(defparameter *default-blocks-per-chunk* (expt 2 16)
  "Number of blocks for the ranks to collectively work on during a single step within APPLY-DISTRIBUTED-GATE.")

(defun apply-distributed-gate (qvm instr &key (blocks-per-chunk *default-blocks-per-chunk*))
  (with-slots (number-of-processes number-of-qubits permutation amplitudes scratch)
      qvm

    (let* ((next-permutation (qubit-permutation instr))
           (addresses (addresses qvm))
           (block-size (block-size addresses))
           (blocks-per-process (blocks-per-process addresses))
           (requests (make-instance 'requests :total (* 2 number-of-processes))))

      (loop :for block-index :from 0 :below (1+ blocks-per-process) :by blocks-per-chunk
            :for start-offset := (* block-index block-size)
            :for end-offset := (* (+ block-index blocks-per-chunk) block-size) :do

              (non-blocking-receive qvm next-permutation start-offset end-offset requests)

              (non-blocking-send qvm next-permutation start-offset end-offset requests)

              ;; XXX Move call to COMPUTE-MATRIX-VECTOR-PRODUCTS here.
              (wait-all requests))

      (compute-matrix-vector-products (quil:gate-matrix instr) scratch amplitudes)

      (update-permutation addresses next-permutation))))

(defun make-arrays (n)
  "Create N empty adjustable arrays."
  (flet ((make-empty-adjustable-array ()
           (make-array 0 :element-type '(unsigned-byte 32) :adjustable t :fill-pointer 0)))

    (make-array n :initial-contents (loop :repeat n :collect (make-empty-adjustable-array)))))

(defun non-blocking-receive (qvm next-permutation start-offset end-offset requests)
  "Iterate over the addresses that should be in the current chunk after applying the next instruction. Start requests to receive the required amplitudes from the ranks where they are currently stored."
  (let* ((addresses (addresses qvm))
         (all-recv-offsets (make-arrays (number-of-processes addresses))))

    (when (< start-offset (number-of-addresses addresses))

      (loop :with next-addresses := (make-addresses-like addresses (rank addresses) next-permutation)
            :for offset :from start-offset :below end-offset :do

              (let ((next-address (get-address-by-offset next-addresses offset)))
                (unless next-address
                  (return))

                (let* ((source-rank (get-rank addresses next-address))
                       (recv-offsets (aref all-recv-offsets source-rank)))
                  (vector-push-extend offset recv-offsets))))

      (post-mpi-irecv qvm all-recv-offsets requests))))

(defun non-blocking-send (qvm next-permutation start-offset end-offset requests)
  "Iterate over all ranks and find which addresses within the current rank are needed by another rank. Aggregate that information, then start sending amplitudes."
  (loop :with addresses := (addresses qvm)
        :with number-of-processes := (number-of-processes addresses)
        :with all-send-offsets := (make-arrays number-of-processes)

        :for count :from 0 :below number-of-processes
        :for target-rank := (mod (+ count (rank qvm)) number-of-processes)
        :for send-offsets := (aref all-send-offsets target-rank)
        :for target-addresses := (make-addresses-like addresses target-rank next-permutation) :do

          (loop :for target-offset :from start-offset :below end-offset :do

            (let ((target-address (get-address-by-offset target-addresses target-offset)))
              (unless target-address
                (return))

              (a:when-let ((source-offset (offset-of addresses target-address)))
                (vector-push-extend source-offset send-offsets))))

        :finally (post-mpi-isend qvm all-send-offsets requests)))
