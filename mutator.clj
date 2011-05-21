(in-ns 'neat.env)
 
;; MUTATION - MUTATOR-FUNCTIONS
;; m: g X gs X e -> gs X e

;; gene mutator functions
(defn param-mutator-f [param-mutator-f]
  (fn [g gs e] 
    (let [ng (assoc g :param (param-mutator-f (:param g)))]
      [(replace {g ng} gs) e])))
(defn add-icon-mutator-f [param-gen-f]	;gene must be node!  --- how I add to inputs?? :O
  (fn [g gs e]
    (let [ns (gs->typ-gs (node-g->rem-connected-node-gs (node-g->other-gs gs g) g) :node)]
      (if (seq ns)
	(let [n (rand-element ns)
	      c (con-gene (:i e) :con (:i n) (:i g) (param-gen-f))
	      e (inc-i e)]
	  [(conj gs c) e])
	[gs e]))))
(defn add-rec-icon-mutator-f [param-gen-f]
  (fn [g gs e]
    (if (not (node-g->rec-connected? gs g))
      (let [c (con-gene (:i e) :con (:i g) (:i g) (param-gen-f))
	    e (inc-i e)]
	[(conj gs c) e])
      [gs e])))
(defn rem-icon-mutator-f []		;gene must be node!
  (fn [g gs e]
    (let [cs (node-g->icon-gs gs g)]
      (if (seq cs)
	(let [c (rand-element (node-g->icon-gs gs g))]
	  [(gs->rem-corresp-g gs c) e])
	[gs e]))))
(defn rem-node-mutator-f []   		;gene must be node!
  (fn [g gs e]
    (let [cs (node-g->iocon-gs gs g)]
      [(gs->rem-corresp-gs gs (conj cs g)) e])))
(defn inject-node-mutator-f [node-param-gen-f con-param-gen-f]
  (fn [g gs e]
    (let [n (node-gene (:i e) :node (node-param-gen-f))
	  e (inc-i e)
	  ic (con-gene (:i e) :con (:ii g) (:i n) (con-param-gen-f))
	  e (inc-i e)
	  oc (con-gene (:i e) :con (:i n) (:oi g) (con-param-gen-f))
	  e (inc-i e)]
      [(conj (gs->rem-corresp-g gs  g) n ic oc) e])))
      

;; COMBINED MUTATOR F
(defn typs-mut-f [typs p mut-f]
  (fn [g gs e]
    (if (some #(= % (:typ g)) typs)
      (if (< (rand 1.0) p)
	(mut-f g gs e)
	[gs e])
      [gs e])))
;; the overall which could be used :)
(defn combined-g-mut-f [typs-p-mut-fs]
  (let [mfs (map (fn [[typs p mf]] (typs-mut-f typs p mf)) typs-p-mut-fs)]
    (fn [g gs e]
      (if (gs->corresp-g gs g)
	(reduce (fn [[ngs ne] mf] (mf g ngs ne)) [gs e] mfs)
	[gs e]))))
(defn combined-gs-mut-f [& typs-p-mut-fs] ;typ p mutf ... format 
  (let [mf (combined-g-mut-f (tupelize typs-p-mut-fs 3))]
    (fn [gs e]
      (reduce (fn [[ngs ne] g] (mf g ngs ne)) [gs e] gs))))

;; combined mutator function
(defn neat-gs-mut-f 
  [input-typs output-typs 
   p-mut-con p-mut-node p-add-con p-add-rec-con p-rem-con p-inject 
   con-param-mut-f con-param-gen-f 
   node-param-mut-f node-param-gen-f]

  (combined-gs-mut-f [:con] p-mut-con (param-mutator-f con-param-mut-f)
		     (conj output-typs :node) p-mut-node (param-mutator-f node-param-mut-f)
		     (conj output-typs :node) p-add-con (add-icon-mutator-f con-param-gen-f)
		     (conj output-typs :node) p-add-rec-con (add-rec-icon-mutator-f con-param-gen-f)
		     (conj output-typs :node) p-rem-con (rem-icon-mutator-f)
;;		     [:node] p-rem-node (rem-node-mutator-f)  ;; ??? use it ??? :)
		     [:con] p-inject (inject-node-mutator-f node-param-gen-f con-param-gen-f))
)

