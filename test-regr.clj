(ns neat.funcregr-test
  (:refer-clojure)
  (:refer msl)
  (:refer msl.seq-numerical)
  (:refer msl.rand)
  (:refer neat.env))


;; the test data
;; (def cases (funcregr/test-sin+cos))
(def cases (funcregr/test-sin))
(def ikeys (range (count (:x (first cases)))))
(def ekeys (map (fn [x] (+ x (count ikeys))) (range (count (:y (first cases))))))
(def template-f (nn.state-funcs/signum--1+1))
(def start-value -1)

(defn vector->ikey-e [x]
  (reduce (fn [h [k v]] (assoc h k v)) 
	  {}
	  (map list ikeys x)))
(defn ekey-e->vector [e]
  (map (fn [k] (get e k)) ekeys))

(defn case->err [c zs0 uf]
  (let [[e zs] (uf (vector->ikey-e (:x c)) zs0)]
    (mse (ekey-e->vector e) (:y c))))
(defn cases->errs [cs zs0 uf]
  (reduce + (map (fn [c] (case->err c zs0 uf)) cs)))       

(defn regr-log-f []
  (fn [gs]
    (if (not (nil? (:id gs))) (println gs))
    (let [[zs0 uf] (nn-genetic/gs->zs0-upd-f gs)]
      (cases->errs cases zs0 uf))))
(defn regr-fit-f []
  (fn [log]
    log))

;; neat enviroment
(defn neat-e []
  (let [[gs i cross-f mut-f dist-f] (nn-genetic/neat-gs-functions 
				     ikeys
				     ekeys 
				     (nn.state-funcs/logistic--1+1 0.1)
				     -1  
				     0.25    
				     0.1  
				     0.1     
				     0.0   
				     0.15  
				     (fn [x] (rand-normal x 0.1)) 
				     (fn [] (rand-normal 0.0 1.0))
				     1.0
				     1.0 
				     3.0 
				     (fn [g0 g1] (cond (nn-genetic/typ? g0 :neuron) (Math/abs (- (:b g0) (:b g1)))
						       (nn-genetic/typ? g1 :con) (Math/abs (- (:w g0) (:w g1))))))						   
	e (neat.env/init-env i  
			     (regr-log-f) 
			     (regr-fit-f)  
			     <  
			     cross-f 
			     mut-f 
			     dist-f 
			     0.75
			     0.6)]
    e))

;; (defn t-evo []
;;   (let [e (neat-e)
;; 	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
;; 	e (neat.env/e-gs->n-mut-add-as-inds e gs 20)]
;;     (e->new-gener-e e))) 
(defn t-evos []
  (let [e (neat-e)
	[gs i] (nn-genetic/ikeys-ekeys->starting-gs (:i e) ikeys ekeys template-f start-value)
	e (neat.env/set-i e i)
	ind (neat.env/gs->ind e gs) 
	[e inds] (neat.env/ind->n-mut-ind e ind 15) 
	inds (neat.env/inds-obj-fit e inds) 
	e (neat.env/add-inds e inds)

	b0 nil
	e0 nil
	b1 nil 
	e1 nil]

    (loop [e e
	   i 0]
      (if (< i 1500)
	(do
	  
	  (let [bi (neat.env/e->best-ind e)]

	    ;; (send b0 (fn [_] @b1))
	    ;; (send b1 (fn [_] bi))
	    ;; (if (:obj-fit @b1) (:obj-fit @b0)

	    (println "best:" (:obj-fit bi) (count (:gs bi)) (count (nn-genetic/gs->neuron-gs (:gs bi))))
	    (recur (e->new-gener-e e) (inc i))))))))
  
 
