;; the normal neat algorithm..
(ns neat
  (:refer-clojure)
  (:refer msl)
  (:refer msl.seq-numerical)
  (:refer msl.rand))


;; the whole function
(defn gs-gs-e->cross-gs [gs-fitter gs e]
  ((:cross-f e) gs-fitter gs))
(defn gs-e->mut-gs-e [gs e]
  ((:mut-f e) gs e))
(defn gs-gs-e->mut-cross-gs-e [gs-fitter gs e]
  (gs-e->mut-gs-e (gs-gs-e->cross-gs gs-fitter gs e) e))


(ns neat.ind
  (:refer msl)
  (:refer nn-genetic)
  (:refer msl.seq-numerical))

;; individual
;; log - data from the individuals performance
;; fit - the fitness of the solution
(defn init-ind [gs]
  {:gs gs
   :obj-log nil
   :obj-fit nil})

(defn obj-log [i obj-log-f]
  (assoc i :obj-log (obj-log-f (:gs i))))
(defn obj-fit [i obj-fit-f]
  (assoc i :obj-fit (obj-fit-f (:obj-log i))))
(defn obj-log-fit [i obj-log-f obj-fit-f]
  (obj-fit (obj-log obj-log-f) obj-fit-f))


;; HANDLING INNOVATION NR
;; new innovation number - iagent - is a counter from the enviroment
(ns neat.env
  (:refer neat)
  (:refer msl)
  (:refer msl.rand)
  (:refer msl.seq-numerical)
  (:refer nn-genetic)
  (:refer clojure.contrib.math))


;; ENVIROMENT
(defn init-env [i obj-log-f obj-fit-f obj-fit-r gs-cross-f gs-mut-f gs-dist-f gs-comp-thr-c specie-succ-rat-c]
  {:i i
   :inds []

   ;; objective function
   :obj-log-f obj-log-f
   :obj-fit-f obj-fit-f
   :obj-fit-r obj-fit-r
   
   ;; mutation etc..
   :gs-cross-f gs-cross-f
   :gs-mut-f gs-mut-f

   ;; speciation functions
   :gs-dist-f gs-dist-f
   :gs-comp-thr-c gs-comp-thr-c
   :specie-succ-rat-c specie-succ-rat-c})

(defn set-i [e i]
  (assoc e :i i))
(defn inc-i [e]
  (assoc e :i (inc (:i e))))

(defn ind-obj-log [e ind]
  (neat.ind/obj-log ind (:obj-log-f e)))
(defn ind-obj-fit [e ind]
  (neat.ind/obj-fit ind (:obj-fit-f e)))
(defn inds-obj-fit [e inds]
  (map #(ind-obj-fit e %) inds))

(defn <ind? [e ind0 ind1]
  ((:obj-fit-r e) (:obj-fit ind0) (:obj-fit ind1)))
(defn inds-fit-sort [e inds]
  (sort #(<ind? e %1 %2) inds))
(defn inds-fit-sum [inds]
  (reduce + (map :obj-fit inds)))
(defn inds-mean-fit [is]
  (avg (map :obj-fit is)))


;; adding new genomes
(defn gs->ind [e gs]
  (ind-obj-log e (neat.ind/init-ind gs)))
(defn gss->ind [e gss]
  (map #(gs->ind e %) gss))
(defn add-ind [e ind]
  (assoc e :inds (conj (:inds e) ind)))
(defn add-inds [e inds]
  (reduce (fn [e ind] (add-ind e ind)) e inds))
(defn gs->add-ind [e gs]
  (add-ind e (gs->ind e gs)))
(defn gss->add-inds [e gss]
  (reduce (fn [e gs] (gs->add-ind e gs)) e gss))
;; (defn add-gs [e gs]
;;   (add-ind e (ind-obj-log e (neat.ind/init-ind gs))))


;; applying crossover and mutation..
(defn gs-gs->cross-gs [e gs0 gs1]
  ((:gs-cross-f e) gs0 gs1))
(defn gs->mut-gs [e gs]
  (let [[gs i] ((:gs-mut-f e) gs (:i e))]
    [(assoc e :i i)
     gs]))

(defn ind->mut-ind [e ind]
  (let [[e gs] (gs->mut-gs e (:gs ind))
	ind (gs->ind e gs)]
    [e ind]))
(defn ind->n-mut-ind [e ind n]
  (reduce (fn [[e is] j] 
	    (let [[e i] (ind->mut-ind e ind)]
	      [e (conj is i)]))
	  [e []]
	  (range n)))
(defn inds->mut-cross-ind [e inds]
  (let [[i0 i1] (inds-fit-sort e (rand-take inds 2))
	gs (gs-gs->cross-gs e (:gs i0) (:gs i1))
	[e gs] (gs->mut-gs e gs)
	i (gs->ind e gs)]
    [e i]))
(defn inds->n-mut-cross-ind [e inds n]
  (reduce (fn [[e is] j]
	    (let [[e i] (inds->mut-cross-ind e inds)]
	      [e (conj is i)]))
	  [e []]
	  (range n)))

;; get the best individual
(defn e-obj-fit-inds [e]
  (assoc e :inds (inds-obj-fit e (:inds e))))
(defn e->best-ind [e]
  (first (inds-fit-sort e (inds-obj-fit e (:inds e)))))
(defn e->best-fit [e]
  (:fit (e->best-ind e)))



;; SPECIATION
;; compatibility distance
;; gs X gs -> n

;; utility functions
(defn ind-ind->comp-dist [e i0 i1]
  ((:gs-dist-f e) (:gs i0) (:gs i1)))
(defn ind-ind->same-specie? [e i0 i1]
  (< (ind-ind->comp-dist e i0 i1) (:gs-comp-thr-c e)))
(defn inds->specie-groups [e inds]
  (group (fn [i0 i1] (ind-ind->same-specie? e i0 i1)) inds))

;; creating the new individual..
(defn specie->new-ind [e inds]
  (let [[e i] (inds->mut-cross-ind e inds)]
    [e i]))
(defn specie->n-new-inds [e inds n]
  (reduce (fn [[e is] j]
	    (let [[e i] (specie->new-ind e inds)]
	      [e (conj is i)]))
	  [e []]
	  (range n)))
(defn specie->complement-to-n [e inds n]
  (let [inds (inds-fit-sort e inds)
	n-inds (count inds)]

    (cond (< n-inds n) (let [[e ninds] (specie->n-new-inds e inds (- n n-inds))]
			 [e (concat inds ninds)])
	  (= n-inds n) [e inds]
	  (< n n-inds) [e (take n inds)])))

;; creating the new specie.. 
(defn specie->new-specie-n [inds avg-fit e]
  (round (/ (inds-fit-sum inds) avg-fit)))
(defn specie->succ [inds e]
  (take (round (* (count inds) (:specie-succ-rat-c e))) (inds-fit-sort e inds)))
(defn specie->succ-n [inds avg-fit e]
  [(specie->succ inds e)
   (specie->new-specie-n inds avg-fit e)])
(defn species->succ-ns [indss avg-fit e]
  (map #(specie->succ-n % avg-fit e) indss))

(defn inc-succ-n [[inds n]]
  [inds (inc n)])
(defn succ-ns->moderate-succ-ns [succ-ns n]
  (let [one-succ-ns (filter (fn [[inds n]] (and (= (count inds) 1) (< 0 n))) succ-ns)
	one-succ-ns (map (fn [[inds n]] [inds 1]) one-succ-ns)
	n-one (count one-succ-ns)
	rem-succ-ns (filter (fn [[inds n]] (< 1 (count inds))) succ-ns)
	n-rem (reduce (fn [i [_ n]] (+ i n)) 0 rem-succ-ns)
	new-rem-succ-ns (loop [succ-ns rem-succ-ns
			       i (- n n-one n-rem)]
			  (if (zero? i)
			    succ-ns
			    (let [ri (rand-int (count succ-ns))]
			      (recur (replace-index succ-ns ri (inc-succ-n (nth succ-ns ri))) (dec i)))))]
    (concat one-succ-ns new-rem-succ-ns)))
	
(defn inds->succ-ns [inds e]
  (let [n (count inds)
	indss (inds->specie-groups e inds)
	avg-fit (inds-mean-fit inds)
	succ-ns (succ-ns->moderate-succ-ns (species->succ-ns indss avg-fit e) n)]
    succ-ns))
(defn succ-n->next-inds [e [inds n]]
  (specie->complement-to-n e inds n))
(defn succ-ns->next-inds [e succ-ns]
  (reduce (fn [[e inds] succ-n]
	    (let [[e s-inds] (succ-n->next-inds e succ-n)]
	      [e (concat inds s-inds)]))
	  [e []] succ-ns))

(defn inds->next-specie-inds [e inds]
  (let [indss (inds->specie-groups e inds)
	succ-ns (inds->succ-ns inds e)
	k (println "species-num:" (count indss))]
    (succ-ns->next-inds e succ-ns)))	

;; the loop
(defn e->new-gener-e [e] 
  (let [e (e-obj-fit-inds e)
	[e inds] (inds->next-specie-inds e (:inds e))]    
    (assoc e :inds inds)))
	

	
 
;; SETUP starting genoms
;; adding the first genom
;; mutate this first genom
;; (defn e->mutate-add-gs-as-ind [e gs]
;;   (let [[gs e] ((:mut-f e) gs e)] 
;;     (e->add-gs-as-ind e gs)))
;; ;; n-times
;; (defn e->mutate-add-gs-as-ind-n [e gs n]
;;   (reduce (fn [e _] (e->mutate-add-gs-as-ind e gs)) e (range n)))


