;;;; qvm-app-ng-tests.asd
;;;;
;;;; Author: appleby

(asdf:defsystem #:qvm-app-ng-tests
  :description "Test Suite for Application server for the QVM."
  :author "Mike Appleby <mappleby@rigetti.com>"
  :license "GNU Affero General Public License v3.0 (See app/LICENSE.txt)"
  :version (:read-file-form "VERSION.txt")
  :depends-on (#:qvm-app-ng
               #:alexandria
               #:drakma
               #:fiasco
               #:lparallel
               #:uiop
               #:yason)
  :perform (asdf:test-op (o s)
                         (uiop:symbol-call ':qvm-app-ng-tests
                                           '#:run-qvm-app-ng-tests))
  :pathname "app-ng/tests/"
  :serial t
  :components ((:file "package")
               (:file "test-validators")
               (:file "test-entry-point")
               (:file "test-rpc-api")
               (:file "test-handlers")
               (:file "test-concurrency")
               (:file "test-safety-hash")
               (:file "test-job")
               (:file "suite")))
