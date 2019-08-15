;;; linear-algebra.lisp
;;;
;;; Author: Juan M. Bello-Rivas

(in-package :dqvm2)

(declaim (optimize (speed 3) (safety 0) (debug 0)))

(defun zgemv (ptr-a m n ptr-x ptr-y)
  "Specialized interface to the ZGEMV BLAS subroutine.

Computes the matrix multiplication of matrix pointed to by PTR-A with the static vector pointed to by PTR-X and stores the result in the static vector pointed to by PTR-Y.

The difference with MAGICL.BLAS-CFFI:%ZGEMV is that this interface allows for pointers to static vectors."
  (cffi:with-foreign-objects ((ptr-m :int32)
                              (ptr-n :int32)
                              (inc :int32)
                              (alpha :double 2)
                              (beta :double 2))
    ;; XXX Make this function return a closure instead, so that these set up steps do not have to be repeated over and over.
    (setf (cffi:mem-ref ptr-m :int32) m
          (cffi:mem-ref ptr-n :int32) n
          (cffi:mem-ref inc :int32) 1
          (cffi:mem-aref alpha :double 0) 1.0d0
          (cffi:mem-aref alpha :double 1) 0.0d0
          (cffi:mem-aref beta :double 0) 0.0d0
          (cffi:mem-aref beta :double 1) 0.0d0)
    ;; Recall that Common Lisp and Fortran use, respectively, row-major and column-major orderings.
    (magicl.blas-cffi::%%zgemv "N" ptr-m ptr-n alpha ptr-a ptr-m ptr-x inc
                               beta ptr-y inc)))

(defun compute-matrix-vector-products (matrix input-array output-array)
  "Loop over blocks and do the matrix-vector multiplications, $y = A x$, where A is MATRIX, INPUT-ARRAY is x, and OUTPUT-ARRAY is y."
  (magicl.cffi-types:with-array-pointers ((ptr-a (magicl::matrix-data matrix)))
    (loop :with stride := (magicl:matrix-cols matrix)
          :with octets-per-cflonum := qvm::+octets-per-cflonum+
          :for i :from 0 :below (length input-array) :by stride
          :do (let* ((offset (* i octets-per-cflonum))
                     (ptr-x (static-vector-pointer input-array :offset offset))
                     (ptr-y (static-vector-pointer output-array :offset offset)))
                (zgemv ptr-a stride stride ptr-x ptr-y)))))
