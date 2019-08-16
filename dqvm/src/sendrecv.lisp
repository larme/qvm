;;; sendrecv.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package #:dqvm2)

(defmethod post-mpi-irecv ((qvm distributed-qvm) all-recv-offsets requests)
  "Post non-blocking MPI_Irecv requests for the amplitudes specified in ALL-RECV-OFFSETS, to be stored in the scratch array of QVM. The associated MPI_Requests are handled by REQUESTS."
  (let ((number-of-processes (number-of-processes qvm))
        (ptr-recvdata (static-vector-pointer (scratch qvm)))
        (tag (qvm::pc qvm)))

    (loop :for source-rank :from 0 :below number-of-processes
          :for recv-offsets := (aref all-recv-offsets source-rank)
          :for count := (length recv-offsets)
          :when (plusp count) :do

            (cffi:with-foreign-array (foreign-recv-offsets recv-offsets `(:array :int ,count))
              ;; XXX store the offsets in a foreign array to begin with.
              (with-mpi-type-indexed-block (recv-type count 1 foreign-recv-offsets +mpi-cflonum+)
                (mpi::%mpi-irecv ptr-recvdata 1 (cffi:mem-ref recv-type 'mpi:mpi-datatype)
                                 source-rank tag +mpi-comm-world+
                                 (get-next-request requests)))))))

(defmethod post-mpi-isend ((qvm distributed-qvm) all-send-offsets requests)
  "Post non-blocking MPI_Isend requests for the amplitudes specified in ALL-SEND-OFFSETS, to be stored in the scratch array of QVM. The associated MPI_Requests are handled by REQUESTS."
  (let ((number-of-processes (number-of-processes qvm))
        (ptr-senddata (static-vector-pointer (amplitudes qvm)))
        (tag (qvm::pc qvm)))

    (loop :for target-rank :from 0 :below number-of-processes
          :for send-offsets := (aref all-send-offsets target-rank)
          :for count := (length send-offsets)
          :when (plusp count) :do

            (cffi:with-foreign-array (foreign-send-offsets send-offsets `(:array :int ,count))
              ;; XXX store the offsets in a foreign array to begin with.
              (with-mpi-type-indexed-block (send-type count 1 foreign-send-offsets +mpi-cflonum+)
                (mpi::%mpi-isend ptr-senddata 1 (cffi:mem-ref send-type 'mpi:mpi-datatype)
                                 target-rank tag +mpi-comm-world+
                                 (get-next-request requests)))))))
