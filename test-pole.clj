(ns neat.test
  (:refer-clojure)
  (:refer msl)
  (:refer msl.seq-numerical)
  (:refer msl.rand)
  (:refer neat.env)
  (:refer neat.genom)
  (:refer neat.repr-nn)
  (:refer pole))
 
(def ityps [:o0 :o1 :o2 :x0 :x1 :x2])
;; (def ityps [:o1 :x1])			;okay.. it's no use =D it's working with every setting.. (but shouldnt :O)
(def otyps [:F])
(defn pole-obj-f []
  (fn [gs] 
    (let [[zs0 uf] (neat.repr-nn/gs->nn-zs0-upd-f gs ityps otyps)]
;;      (println "sim" (pole/simulate-controller zs0 uf 2.5))
      (pole/simulate-controller zs0 uf 5.0))))

;; !!!!!
;; PROBLEM!!! :O -> i have the inputs  (like om0) will be refreshed by net-update  and exported to the env :(((
;; !!!!! 
;; THIS SUCKS ELEMENTARY =D

;; polebalance - om0 om1 om2 x0 x1 x2 
(defn neat-e []
  (neat.env/init-env  
   []  
   0
   (pole-obj-f)
   >  
   (neat.genom/crossover-f 0.5)
   (neat.env/neat-gs-mut-f ityps
			   otyps
			   0.2
			   0.1  
			   0.05    
			   0.0 		;rec con  
			   0.025    
			   0.025  
			   (fn [x] (rand-normal x (Math/abs (/ x 50))))
			   (fn [] (- 1 (rand 2.0))) 
			   (fn [x] (rand-normal x (Math/abs (/ x 50))))
			   (fn [] (- 1 (rand 2.0))))
   (neat.env/compat-dist-f 1.0    
			   1.0  
			   3.0 
			   (fn [x0 x1] (Math/abs (- x0 x1)))) ;maybe bigger param diff function?? 
   4.0 
   0.6))
       

(defn t-evo []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	e (neat.env/e-gs->n-mut-add-as-inds e gs 20)]
    (e->new-gener-e e)))
(defn t-evos []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	e (neat.env/e-gs->n-mut-add-as-inds e gs 10)]
    (loop [e e
	   i 0]
      (if (< i 50)
	(do
	  (println "best" (:id (:best-ind e)) (:fit (:best-ind e)) (count (:gs (:best-ind e)) ))
	  (recur (e->new-gener-e e) (inc i)))))))
  


;;  (def e1 (neat.repr-nn/e->add-starting-gss e0 10 ityps otyps))
     
(defn t []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	[zs uf] (neat.repr-nn/gs->nn-zs0-upd-f gs ityps otyps)]
    (uf (pole/init-e) zs)))
(defn t2 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	of (pole-obj-f)]
    (of gs)))
(defn t3 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	mf (neat.env/neat-gs-mut-f ityps
				   otyps
				   1 1 1 1 1 1 
				   ;; 0.2 
				   ;; 0.1  
				   ;; 0.05 
				   ;; 0.0 		;rec con  
				   ;; 0.025  
				   ;; 0.025 
				   (fn [x] (rand-normal x (Math/abs (/ x 50)))) 
				   (fn [] (- 1 (rand 2.0))) 
				   (fn [x] (rand-normal x (Math/abs (/ x 50))))
				   (fn [] (- 1 (rand 2.0))))]
    (mf gs e)))
(defn t4 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	[e i] (neat.env/e-gs->mut-add-as-ind e gs)]
    e))
(defn t5 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	e (neat.env/e-gs->n-mut-add-as-inds e gs 5)]
    e))
(defn t6 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	cong (nth gs 10)
	nodg (nth gs 7)
	
	pmf (fn [x] (rand-normal x (Math/abs (/ x 50))))
	pgf (fn [] (- 1 (rand 2.0))) 
	
	m0 (neat.env/param-mutator-f pmf)
	m1 (neat.env/rem-icon-mutator-f)
	m2 (neat.env/add-icon-mutator-f pgf)
	m3 (neat.env/add-rec-icon-mutator-f pgf)
	m4 (neat.env/inject-node-mutator-f pgf pgf)
	m5 (neat.env/rem-node-mutator-f)]
   
    (println 
     ;; (m0 cong gs e)
     ;; (m0 nodg gs e)
     ;; (m1 nodg gs e)
     ;; (m2 nodg gs e)
     ;; (m3 nodg gs e)
     ;; (m4 cong gs e)
     ;; (m5 nodg gs e)
)))
(defn t7 []
  (let [e (neat-e)
	[e gs] (neat.repr-nn/starting-gs e ityps otyps (fn [] (- 1 (rand 2.0))))
	e (neat.env/e-gs->n-mut-add-as-inds e gs 10)]

;;    (in-specie? (first (:inds e)) (second (:inds e)) e)
;;    (println (map :param (map first (map :gs (:inds e)))))
    (map (fn [i0] (ind-ind->comp-dist i0 (first (:inds e)) e)) (rest (:inds e)))
    ;; (doseq [i0 (rest (:inds e))]
    ;;   (let [k (println "jijj")
    ;; 	    d (ind-ind->comp-dist i0 (first (:inds e)) e)]
    ;; 	(println "aijj" d))) 
;;    (println "AAA")
    ;; ITS POSSIBLE THAT THE PROBLEM IS AROUND THE mutation!!! :O
))


    ;; (in-specie? (first (:inds e)) (second (:inds e)) e)))


 