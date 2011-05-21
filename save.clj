;; the normal neat algorithm..
(ns neat
  (:refer-clojure)
  (:refer msl)
  (:refer msl.seq-numerical)
  (:refer msl.rand))

;; a few common gene typs..!!! 
(ns neat.gene
  (:refer neat))
(defn base [i typ]
  {:i i 				;innovation numbers
   :typ typ})

;; connection gene
(defn con-gene [i iid oid w]
  (assoc (base i :con) 
    :iid iid
    :oid oid
    :w w))

;; ;; neat starting env
;; (defn mk-starting-env [starting-pop 

;; genom
(ns neat.genom
  (:refer neat)
  (:refer msl)
  (:refer msl.seq-numerical)
  (:refer msl.rand))

;; sorting by it's innovation number
(defn innov-sort [gs]
  (sort (fn [g0 g1] (< (:i g0) (:i g1))) gs))

;; lineup - matching genes
(defn common-disjoint-excess [gs0 gs1]
  (loop [gs [gs0 gs0]
	 cs []
	 ds []
	 es []]
    (let [[[a0 & ags0] [a1 & ags1]] gs]
      (cond (not (and a0 a1)) [cs ds es]	
	    (or (nil? a0) (nil? a1)) (recur [ags0 ags1] cs ds (conj es [a0 a1]))
	    (= (:i a0) (:i a1)) (recur [ags0 ags1] (conj cs [a0 a1]) ds es)
	    (< (:i a0) (:i a1)) (recur [ags0 (cons a1 ags1)] cs (conj ds [a0 nil]) es)
      	    (> (:i a0) (:i a1)) (recur [(cons a0 ags0) ags1] cs (conj ds [nil a1]) es)))))


;; crossing over function
(defn crossover-f [prob-choose-fitter] ;0.5 in normal crossover
  (fn [gs-fitter gs]
    (let [[cs ds es] (common-disjoint-excess gs-fitter gs)]
      (innov-sort (concat (map (fn [[g0 g1]] (if (< (rand) prob-choose-fitter) g0 g1)) cs)
			  (remove nil? (map first ds))
			  (remove nil? (map first es)))))))
	      


;; mutation    
;; all gene-param-mutate-f ->  
(defn param-mutate-f [gene-param-mutate-fs]
(defn if-typ-f [typ f]
  (fn [g] (if (= (:typ g) typ)
	    (f g))))

;; distance evolution we should introduce?? :O
(defn gene-mutate-fixed-var-f [prob-mutate-gene field variance]
  (fn [g] (if (< (rand) prob-mutate-gene)
	    (assoc g field (rand-normal (field g) variance))
	    g)))

;; it's applying all -> rule should apply themselves only to their typs
(defn gene-rules-mutate-f [mutate-fs]
  (fn [g] 
    (reduce (fn [f] (f g)) mutate-fs)))
(defn genom-rules-mutate-f [mutate-f]
  (fn [gs]
    (map (fn [g] (mutate-f g)) gs)))
	 

;; speciation	  
;; compatibility distance
(defn compatibility-distance [gs0 gs1 c-disjoint c-excess c-param-diff param-diff-f]
  (let [[cs ds es] (common-disjoint-excess gs0 gs1)
	[n d e] (map count [cs ds es])
	avg-diff (avg (map (fn [c0 c1] (param-diff-f c0 c1)) cs))]
    
    (+ (* c-disjoint (/ d n))
       (* c-excess (/ e n))
       (* c-param-diff avg-diff))))
    
(defn speciation [gs comp-dis-f comp-dis-thres]
  (group (fn [g0 g1] (< (comp-dis-f g0 g1) comp-dis-thres))
	 gs))




(ns neat.individual
  (:refer neat))

(defn init [genes]
  {:gs genes
   :obj nil				; objective-score
   :sh-obj nil})			; shared-objective-score



(ns neat.env
  (:refer neat)
  (:refer msl)
  (:refer msl.rand)
  (:refer clojure.contrib.math))

;; neat individual

;; two types of mutate
;; gene - their parameters change - no number increment
;; structural - new innovation number to the parts 

(defn init 
  [individums 
   objective-f 
   gene-mutation-f 
   structural-mutation-f

   const-disjoint
   const-excess
   const-paramdiff

   const-compatibility-threshold
   const-success-ratio 			;??
   ]

  {:inds individums

   :obj-f objective-f
   :g-mut-f gene-mutation-f		; g -> g
   :s-mut-f structural-mutatation-f	; g -> [g g ..]

   :cppm const-prob-param-mutation
   :cpsm const-prom-structural-mutation

   :cd const-disjoint
   :ce const-excess
   :cpd const-param-distance
   :pdf param-distance-f
   
   :cct const-compatibility-threshold
   :csr const-success-ratio 

   :i 0}
)

(defn compatibility-dist [i0 i1 e]
  (neat.genom/compatibility-distance (:gs i0) (:gs i1) (:cd e) (:ce e) (:cpd e) (:pdf e)))
(defn speciation [e])

		  
;; individual

(defn init-working-env []
  {:inds individuals
   
   :obj-f objective-f			; gs -> R    (= fitness :))
   :obj-cmp-f objective-comparator-f	; f X f -> bool

   :offs-f offspring-f			; gs X gs X env -> gs
   :comp-f compatibility-f		; gs X gs -> bool
   
   :i 0})				; innovation number

;; group the individuals into species
(defn speciate-inds [is comp-f]   
  (group (fn [i0 i1] (comp-f (:gs i0) (:gs i1))) is))

;; give the objective function scores for the individuals
(defn obj-inds [is obj-f]
  (map (fn [i] (assoc i :obj (obj-f (:gs i)))) is))
(defn sh-obj-inds-in-specie [is]
  (let [n (count is)]
    (map (fn [i] (assoc i :sh-obj (/ (:obj i) n))) is)))
(defn sh-obj-inds-in-species [speciated-is]
  (map (fn [is] (sh-obj-inds-in-group is)) speciated-is))

;; sorting for mating
(defn obj-sort-specie [is obj-cmp-f]
  (sort (fn [i0 i1] (obj-cmp-f (:obj i0) (:obj i1))) is))
(defn obj-sort-species [speciated-is obj-cmp-f]
  (map (fn [is] (obj-sort-specie obj-cmp-f)) speciated-is))
;;  (sort (fn [[i0 & is0] [i1 & is1]] (obj-cmp-f (:obj i0) (:obj i1))) speciated-is))

;; select for mating
(defn select-specie-fittest [is select-rat]
  (let [n (count is)]
    (if (= n 1)
      [n is]
      (let [nsel0 (ceil (* n select-rat))
	    nsel (if (< nsel0 2) 2 nsel0)]
	[n (take nsel is)]))))
(defn select-species-fittest [speciated-is select-rat]
  (map (fn [is] (select-species-fittest is select-rat))))

;; mating of individuals
(defn next-specie-genom-generation [n is offspring-f]
  (loop [j (- n (count is))
	 ngs []]
    (if (zero? j)
      ngs
      
      (let [[i0 i1] (rand-take is 2)]
	(recur (dec j) (conj (offspring-f i0 i1))))))) ;problem here because of the innovation nr :O
(defn next-specie-generation [n is offspring-f obj-f]
  (concat is (for [gs (next-specie-genom-generation n is offspring-f)]
	       {:gs gs
		:obj (obj-f gs)})))
		
	      
