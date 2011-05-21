
(ns neat.repr-nn
  (:refer-clojure)
  (:refer msl)
  (:refer msl.rand)
  (:refer neat.gene)
  (:refer neat.genom)
  (:refer neat.env))

;; STARTING-SETUP gene representation
(defn typ->ng [e typ w-rand-f]
  (let [ng (node-gene (:i e) typ (w-rand-f))]
    [(inc-i e) ng]))
(defn typs->ngs [e typs w-rand-f]
  (reduce (fn [[e ngs] t] 
	    (let [[e ng] (typ->ng e t w-rand-f)]
	      [e (conj ngs ng)]))
	  [e []]
	  typs))
(defn cg-between [e ii oi w-rand-f]
  [(inc-i e) (con-gene (:i e) :con ii oi (w-rand-f))])
(defn cg-between-ngs [e ing ong w-rand-f]
  (cg-between e (:i ing) (:i ong) w-rand-f))
(defn cgs-between-ngs [e ings ongs w-rand-f] ;every possible pairing
  (reduce (fn [[e cgs] [ing ong]]
	    (let [[e cg] (cg-between-ngs e ing ong w-rand-f)]
	      [e (conj cgs cg)]))
	  [e []]
	  (for [ing ings ong ongs] (vector ing ong))))

;; starting gs
;; every output will receive a fresh node with inputs etc..
(defn otyp->ngs-cgs [e otyp ings w-rand-f]
  (let [[e ong] (typ->ng e otyp w-rand-f)
	[e ng] (typ->ng e :node w-rand-f)
	[e icgs] (cgs-between-ngs e ings [ng] w-rand-f)
	[e ocg] (cg-between-ngs e ng ong w-rand-f)]
    [e (conj icgs ng ocg ong)]))
(defn otyps->ngs-cgs [e otyps ings w-rand-f]
  (reduce (fn [[e gs] otyp] 
	    (let [[ne ngs] (otyp->ngs-cgs e otyp ings w-rand-f)]
	      [ne (concat gs ngs)]))
	  [e []]
	  otyps))

(defn rand1 []
  (- 1 (rand 2)))
(defn starting-gs [e input-typs output-typs w-rand-f]
  (let [[e ings] (typs->ngs e input-typs w-rand-f)
	[e gs] (otyps->ngs-cgs e output-typs ings w-rand-f)]
    [e (concat gs ings)]))
	


;; translate GENES to INFS CINFS
;; index hash - gs-i -> node-index(in the state vector)
(defn ngs->index-hash [ngs] 
  (reduce (fn [h [ng i]] (assoc h (:i ng) i))
	  {}
	  (for [ng ngs 
		i (range (count ngs))] 
	    (vector ng i))))
(defn cg->cinf [cg index-hash]
  (nn/con-inf (get index-hash (:oi cg)) (get index-hash (:ii cg)) (:param cg)))
(defn ng->ninf [ng index-hash]
  (nn/node-inf (get index-hash (:i ng)) (:typ ng) (:param ng)))
(defn gs->ninfs-cinfs [gs]
  (let [ngs (gs->node-gs gs)
	cgs (gs->con-gs gs)
	ih (ngs->index-hash ngs)
	ninfs (map (fn [ng] (ng->ninf ng ih)) ngs)
	cinfs (map (fn [cg] (cg->cinf cg ih)) cgs)]
    [ninfs cinfs]))
(defn gs->nn-zs0-upd-f [gs input-typs output-typs zf]
  (let [[ninfs cinfs] (gs->ninfs-cinfs gs)]
    (nn/ninfs-cinfs->zs0-upd-f ninfs cinfs input-typs output-typs zf)))

	
;; STARTING NEAT ENV
(defn e->add-starting-gss [e n input-typs output-typs]
  (let [[e gs] (starting-gs e input-typs output-typs rand1)]
    (e-gs->n-mut-add-as-inds e gs n)))




