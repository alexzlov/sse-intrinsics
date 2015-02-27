(defpackage "MANDELBROT"
  (:use "CL"
        "SB-EXT"          ;; for SIMD-PACK
        "SB-C"))          ;; for DEFKNOWN etc

(in-package "MANDELBROT")

(defknown (f4+ f4* f4-) ((simd-pack single-float) (simd-pack single-float))
  (simd-pack single-float)
  (movable flushable always-translatable)
  :overwrite-fndb-silently t)

(defknown f4<= ((simd-pack single-float) (simd-pack single-float))
    (simd-pack (signed-byte 32))
    (movable flushable always-translatable)
  :overwrite-fndb-silently t)

(defknown i4- ((simd-pack (signed-byte 32)) (simd-pack (signed-byte 32)))
    (simd-pack (signed-byte 32))
    (movable flushable always-translatable)
  :overwrite-fndb-silently t)

(defknown f4-sign-mask (simd-pack) (unsigned-byte 4)
    (movable flushable always-translatable)
  :overwrite-fndb-silently t)


(in-package "SB-VM")

;; new VOP definition.
;; when there is a call to mandelbrot::f4+ that can be used
;; both in fast and in safe code.
(define-vop (mandelbrot::f4+)
  (:translate mandelbrot::f4+)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (single-sse-reg)))
  (:result-types simd-pack-single)
  (:generator 4
              (cond ((location= r y)
                     (inst addps y x))
                    (t
                     (move r x)
                     (inst addps r y)))))


(define-vop (mandelbrot::f4*)
  (:translate mandelbrot::f4*)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (single-sse-reg)))
  (:result-types simd-pack-single)
  (:generator 4
              (cond ((location= r y)
                     (inst mulps y x))
                    (t
                     (move r x)
                     (inst mulps r y)))))


(define-vop (mandelbrot::f4-)
  (:translate mandelbrot::f4-)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (single-sse-reg) :from (:argument 0)))
  (:result-types simd-pack-single)
  (:generator 4
              (move r x)
              (inst subps r y)))


(define-vop (mandelbrot::f4<=)
  (:translate mandelbrot::f4<=)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (int-sse-reg) :from (:argument 0)))
  (:result-types simd-pack-int)
  (:generator 4
              (move r x)
              (inst cmpps :le r y)))


(define-vop (mandelbrot::f4-sign-mask)
  (:translate mandelbrot::f4-sign-mask)
  (:policy :fast-safe)
  (:args (x :scs (int-sse-reg single-sse-reg double-sse-reg)))
  (:arg-types simd-pack)
  (:results (r :scs (unsigned-reg)))
  (:result-types unsigned-num)
  (:generator 4
              (inst movmskps r x)))


(define-vop (mandelbrot::f4<=)
  (:translate mandelbrot::f4<=)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (int-sse-reg) :from (:argument 0)))
  (:result-types simd-pack-int)
  (:generator 4
              (move r x)
              (inst cmpps :le r y)))


(define-vop (mandelbrot::i4-)
  (:translate mandelbrot::i4-)
  (:policy :fast-safe)
  (:args (x :scs (int-sse-reg) :target r) (y :scs (int-sse-reg)))
  (:arg-types simd-pack-int simd-pack-int)
  (:results (r :scs (int-sse-reg) :from (:argument 0)))
  (:result-types simd-pack-int)
  (:generator 4
              (move r x)
              (inst psubd r y)))


(define-vop (mandelbrot::f4-sign-mask)
  (:translate mandelbrot::f4-sign-mask)
  (:policy :fast-safe)
  (:args (x :scs (int-sse-reg single-sse-reg double-sse-reg)))
  (:arg-types simd-pack)
  (:results (r :scs (unsigned-reg)))
  (:result-types unsigned-num)
  (:generator 4
              (inst movmskps r x)))


(in-package "MANDELBROT")

(macrolet ((define-stub (name)
             `(defun ,name (x y)
                (,name x y))))
  (define-stub f4+)
  (define-stub f4-)
  (define-stub f4*)
  (define-stub f4<=)
  (define-stub i4-))

(defun f4-sign-mask (x)
  (f4-sign-mask x))


(deftype f4 ()
  '(simd-pack single-float))

(declaim (inline %mandelbrot-iter replicate-float))

(defun replicate-float (x)
  (%make-simd-pack-single x x x x))

(defun %mandelbrot-iter (zr zi cr ci)
  (declare (optimize speed)
           (type f4 zr zi cr ci))
  (let* ((r (f4- (f4* zr zr)
                 (f4* zi zi)))
         (i/2 (f4* zr zi))
         (i   (f4+ i/2 i/2)))
    (values (f4+ r cr)
            (f4+ i ci))))

(declaim (inline %norm^2))

(defun %norm^2 (r i)
  (declare (optimize speed) (type f4 r i))
  (f4+ (f4* r r) (f4* i i)))
