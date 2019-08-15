;;; addresses.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package :dqvm2)

;;; What do we need?

(defvar *print-addresses* t
  "Print address numbers when printing objects of class ADDRESSES.")

(defparameter *default-block-size* 4
  "Default block size for address tables.")

(defclass addresses ()
  ((number-of-qubits          ; XXX Find a proper way to share this with QVM.
    :reader number-of-qubits
    :initarg :number-of-qubits
    :type qvm::non-negative-fixnum
    :initform (error-missing-initform :number-of-qubits)
    :documentation "Number of qubits in the QVM.")
   (number-of-addresses
    :reader number-of-addresses
    :type qvm::non-negative-fixnum
    :documentation "Length of the subvector of the wavefunction indexed by this table")
   (total-number-of-addresses
    :reader total-number-of-addresses
    :type qvm::non-negative-fixnum
    :documentation "Length of the wavefunction")
   (rank                                ; XXX Ditto.
    :reader rank
    :initarg :rank
    :initform (mpi-comm-rank)
    :documentation "MPI rank in charge of this instance.")
   (number-of-processes                 ; XXX Ditto.
    :accessor number-of-processes
    :initarg :number-of-processes
    :initform (mpi-comm-size)
    :documentation "Number of parallel MPI ranks.")
   (block-size
    :reader block-size
    :initarg :block-size
    :type qvm::non-negative-fixnum
    :initform *default-block-size*
    :documentation "Minimum number of amplitudes per single block matrix-vector multiplication. This must be an even number and is typically chosen to be equal to two times the arity of the instruction with the maximum arity in the program to be executed.")
   (blocks-per-process
    :reader blocks-per-process
    :type qvm::non-negative-fixnum
    :documentation "Regular blocks per process (i.e., without including a possibly remaining block).")
   (remainder-blocks
    :reader remainder-blocks
    :type qvm::non-negative-fixnum
    :documentation "Total number of remainder blocks.")
   (total-blocks
    :reader get-total-blocks            ; XXX remove the get- prefix
    :type qvm::non-negative-fixnum
    :documentation "Total number of blocks among all address tables.")
   (permutation
    :accessor permutation
    :initarg :permutation
    :initform nil
    :type permutation
    :documentation "Last qubit permutation evaluated, stored in a format suitable for use by APPLY-QUBIT-PERMUTATION."))

  (:documentation "Table of addresses handled by a single rank of a distributed QVM. The addresses within the table are computed on demand, rather than being explicitly stored."))

(defmethod initialize-instance :after ((addresses addresses) &key)
  (with-slots (number-of-qubits
               total-number-of-addresses number-of-addresses
               rank number-of-processes
               block-size blocks-per-process remainder-blocks
               total-blocks)
      addresses
    (assert (not (minusp number-of-qubits)))
    (assert (not (minusp rank)))
    (assert (plusp number-of-processes))
    (assert (< rank number-of-processes))
    (assert (zerop (mod block-size 2)))

    (multiple-value-bind (bpp rem) (%blocks-per-process addresses)
      (setf total-number-of-addresses (ash 1 number-of-qubits)
            number-of-addresses (%number-of-addresses addresses)
            total-blocks (%get-total-blocks addresses)
            blocks-per-process bpp
            remainder-blocks rem))))

(defun make-addresses-like (addresses new-rank &optional (new-permutation nil new-permutation-p))
  "Create new address table like ADDRESSES. Set the rank (and optionally the permutation) of the new instance to be NEW-RANK (and, respectively, NEW-PERMUTATION)."
  (with-slots (number-of-processes number-of-qubits block-size permutation)
      addresses
    (make-instance 'addresses :rank new-rank
                              :number-of-processes number-of-processes
                              :number-of-qubits number-of-qubits
                              :block-size block-size
                              :permutation (if new-permutation-p
                                               new-permutation
                                               (permutation addresses)))))

(qvm::defun-inlinable %get-total-blocks (addresses)
  "Return the total number of block matrices to consider in the simulation."
  (let ((total-amplitudes (ash 1 (number-of-qubits addresses))
                          ;; Call NUMBER-OF-ADDRESSES instead?
                          ))
    (/ total-amplitudes (block-size addresses))))

(defmethod %blocks-per-process ((addresses addresses))
  "Return the number of regular blocks per process in ADDRESSES and the number of remaining blocks.

Note that this may not be the number of blocks in the processes. To obtain that, you must check the rank and the remaining blocks."
  (floor (%get-total-blocks addresses)
         (number-of-processes addresses)))

(defmethod get-initial-address ((addresses addresses) address)
  "Apply the inverse of the latest permutation in ADDRESSES to ADDRESS."
  (apply-inverse-qubit-permutation (permutation addresses) address))

(defmethod get-rank ((addresses addresses) address)
  "Find the MPI rank of ADDRESS using the table ADDRESSES.

This method returns the rank where the address is located regardless of whether that rank coincides with that of ADDRESSES. If ADDRESS is NIL, the return value is also NIL."
  (when address
    (with-slots (number-of-processes block-size blocks-per-process remainder-blocks total-blocks)
        addresses

      (assert (zerop (mod block-size 2)))

      (let ((initial-address (get-initial-address addresses address)))
        (if (< initial-address (* (- total-blocks remainder-blocks) block-size))
            (nth-value 0 (floor initial-address (* blocks-per-process block-size)))
            (mod (floor initial-address block-size) number-of-processes))))))

(defmethod member-p ((addresses addresses) address)
  "Return T if ADDRESS is in the table ADDRESSES."
  (= (get-rank addresses address) (rank addresses)))

(defmethod offset-of ((addresses addresses) address)
  "Find the index (offset) of the amplitude corresponding to ADDRESS. Return NIL if address is not in ADDRESSES."
  ;; XXX blend member-p with this function to avoid redundant work.
  (when (member-p addresses address)
    (with-slots (block-size blocks-per-process remainder-blocks total-blocks permutation rank)
        addresses

      (let* ((initial-address (get-initial-address addresses address))
             (start-of-remainder-addresses (* (- total-blocks remainder-blocks) block-size)))

        (if (< initial-address start-of-remainder-addresses)
            (- initial-address (* rank blocks-per-process block-size))
            (+ (* blocks-per-process block-size)
               (- initial-address start-of-remainder-addresses (* rank block-size))))))))

(defmethod %number-of-addresses ((addresses addresses))
  "Number of addresses handled by the table ADDRESSES."
  (multiple-value-bind (blocks-per-process remainder)
      (floor (%get-total-blocks addresses)
             (number-of-processes addresses))

    (* (block-size addresses)
       (+ blocks-per-process (if (< (rank addresses) remainder)
                                 1
                                 0)))))

(defmethod effective-permutation ((addresses addresses) next-permutation)
  "Return permutation required to map the original ADDRESSES to NEXT-PERMUTATION.

If we have a table that is the image of (0, 1, ...) by the permutation πₖ and we want to map the original indices to πₕ, then we need to apply the effective permutation πₕ ∘ πₖ⁻¹."
  (let ((curr-permutation (permutation addresses)))
    (append (reverse curr-permutation) next-permutation)))

(defmethod get-effective-permutation ((addresses addresses) next-permutation)
  "Get the next permutation needed to apply PERMUTATION on ADDRESSES.

Returns the effective permutation (i.e., π₂ ∘ π₁⁻¹) as well as the next permutation (i.e., π₂). These permutations are in a format ready to be passed to APPLY-QUBIT-PERMUTATION."
  (let ((prev-permutation (permutation addresses)))
    (append (reverse prev-permutation) next-permutation)))

(defmethod update-permutation ((addresses addresses) next-permutation)
  "Update permutation in ADDRESSES by NEXT-PERMUTATION."
  (setf (permutation addresses) next-permutation)
  addresses)

(defun %get-address-ranges (addresses)
  (with-slots (rank number-of-processes block-size blocks-per-process remainder-blocks)
      addresses

    (let ((total-blocks (get-total-blocks addresses)))

      ;; (format t "blocks-per-process ~a~tremainder ~a~tblock-size ~a~%"
      ;;         blocks-per-process remainder block-size)

      (let* ((k (* block-size blocks-per-process))
             (result (list (* rank k) (* (1+ rank) k))))
        (when (< rank remainder-blocks)
          (let ((n-times-k (* number-of-processes k)))
            (nconc result (list (+ n-times-k (* rank block-size))
                                (+ n-times-k (* (1+ rank) block-size))))))
        (values-list result)))))

(defmethod get-address-by-offset ((addresses addresses) offset)
  "Return the address located at OFFSET within ADDRESSES or NIL if the offset exceeds the number of addresses held."
  (check-type offset (unsigned-byte 64))

  (with-slots (rank number-of-addresses number-of-processes
               block-size blocks-per-process permutation)
      addresses

    (assert (not (minusp offset)))

    (when (< offset number-of-addresses)
      (let ((base-address (if (< offset (* blocks-per-process block-size))
                              (* rank blocks-per-process block-size)
                              (* (+ rank
                                    (* (1- number-of-processes) blocks-per-process))
                                 block-size))))
        (apply-qubit-permutation permutation (+ offset base-address))))))

(defmacro do-addresses ((var addresses &optional result) &body body)
  "Iterate over ADDRESSES using VAR.

The addresses are generated on the fly based on the rank of ADDRESSES and the permutation."
  (a:once-only (addresses)
    (a:with-gensyms (start-0 stop-0 start-1 stop-1 permutation index)
      `(multiple-value-bind (,start-0 ,stop-0 ,start-1 ,stop-1)
           (%get-address-ranges ,addresses)

         (let ((,permutation (permutation ,addresses)))

           (loop :for ,index :from ,start-0 :below ,stop-0
                 :for ,var := (apply-qubit-permutation ,permutation ,index)
                 :do ,@body)

           (when (and ,start-1 ,stop-1)
             (loop :for ,index :from ,start-1 :below ,stop-1
                   :for ,var := (apply-qubit-permutation ,permutation ,index)
                   :do ,@body))

           ,result)))))

(defmethod print-object ((addresses addresses) stream)
  (let ((*print-readably* nil)
        (*print-pretty* nil))
    (print-unreadable-object (addresses stream :type t :identity t)
      (with-slots (rank number-of-processes number-of-qubits
                   block-size permutation)
          addresses

        (format stream "~a ~d ~a ~d ~a ~d ~a ~d ~a ~a"
                (prin1-to-string :rank) rank
                (prin1-to-string :number-of-processes) number-of-processes
                (prin1-to-string :number-of-qubits) number-of-qubits
                (prin1-to-string :block-size) block-size
                (prin1-to-string :permutation) permutation)

        (when *print-addresses*
          (let (address-list)
            (do-addresses (address addresses)
              (push address address-list))
            (when address-list
              (format stream " ~a (~{~d~^ ~})"
                      (prin1-to-string :addresses)
                      (nreverse address-list)))))))))

(defmethod find-next-addresses ((addresses addresses) next-permutation)
  ;; XXX Review docstring.
  "Find amplitude addresses to exchange when applying NEXT-PERMUTATION and the rank where the amplitudes are located.

Returns two vectors: one contains the new addresses and the other contains ranks where said addresses are located.

Example
-------

Consider the following example program:

    I 1
    I 2

There are 3 qubits (0, 1, and 2) and 2^3 = 8 addresses and amplitudes. If we have 2 ranks available, then each of them keeps 2³ / 2 = 4 address/amplitude pairs and the the data is split among them as shown below:

      rank #0 | rank #1
      --------+--------
      0 1 2 3 | 4 5 6 7

Suppose we have already executed \"I 1\" and we are about to execute \"I 2\". Moreover, imagine we are sitting at rank #0, then we must find out two pieces of information:

  1. Which address/amplitude pairs are going to be sent and where are they going (i.e., which rank and which place within the receiving rank).

  2. Which addresses and amplitudes are we going to be receiving from where.

In order to find out the previous pieces of information, it is convenient to keep the following diagram in mind:

      rank #0 | rank #1
      --------+--------
  π₀  0 1 2 3 | 4 5 6 7
  π₁  0 2 1 3 | 4 6 5 7  <--- Addresses after executing \"I 1\"
  π₂  0 4 2 6 | 1 5 3 7  <--- Addresses after executing \"I 2\"

where πₖ is the permutation associated with the k-th instruction and π₀ is the identity.

To see where addresses 0, 2, 1, and 3 are migrating to, we compute

  π₂⁻¹(0, 2, 1, 3) = (0, 2, 4, 6) = (0, 0, 1, 1) mod 2³/2

which tells us that addresses/amplitudes 0 and 2 remain in rank #0 (just take the remainder of the division by 2³/2 to see it) and, similarly, that 1 and 3 are to be sent to rank #1. Note that we also know the order into rank #1's address space (addresses 1 and 3 appear in the 4th and 6th global address or,equivalently, offsets 0 and 2 of rank #1's address space).

Now it remains to find out the addresses that rank #0 will be receiving. This follows by applying

  π₂ ∘ π₁⁻¹(0, 2, 1, 3) = (0, 4, 2, 6),

and then

  π₁⁻¹(0, 4, 2, 6) = (0, 4, 1, 5) = (0, 1, 0, 1) mod 2³/2

tells us that addresses 4 and 6 will come from rank #1 (and in which order)."
  (check-type next-permutation permutation)

  (let ((permutation (permutation addresses))
        (effective-permutation (get-effective-permutation addresses next-permutation))
        (orig-addresses nil)
        (next-addresses nil)
        (target-addresses nil)
        (source-addresses nil))

    (do-addresses (address addresses (values (nreverse orig-addresses)
                                             (nreverse next-addresses)
                                             (nreverse target-addresses)
                                             (nreverse source-addresses)))
      ;; First, find the next address $y = π₂ ∘ π₁⁻¹(x)$, where $x$ is the current address. Next, compute $z = π₁^{-1}(y)$.
      (let* ((next-address (apply-qubit-permutation effective-permutation address))
             (target-address (apply-inverse-qubit-permutation next-permutation address))
             (source-address (apply-inverse-qubit-permutation permutation next-address)))

        (push address orig-addresses)
        (push next-address next-addresses)
        (push target-address target-addresses)
        (push source-address source-addresses)))))
