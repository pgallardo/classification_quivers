clearAll
needsPackage ("ThinSincereQuivers")
needsPackage("NormalToricVarieties")
needsPackage("Polyhedra")
needsPackage("Graphs")
needsPackage("LatticePolytopes")

-- Define the quiver associated to the maximal quiver that generates
-- all toric quivers of dimension equal to two.

f = "ClassificationQuiver2D" << ""

Q2D = {{0,2},{0,3},{0,4},{1,2},{1,3},{1,4}}
K  = toricQuiver(Q2D)
CG = coneSystem(K)
RT = referenceThetas CG

f << length RT << endl

f << RT << endl
--RT = {{-5, -1, 2, 2, 2}}

-- Initialize a list to hold all IsoCases for each i

AllIsoCases = new MutableList from apply(length RT, i -> {});

for i from 0 to (length RT - 1) do {
    wa = RT#i;
    fa = incInverse(wa, K) + {-7, -11, 18, 7, 11, -18};
    QuiverA = toricQuiver(Q2D, fa);
    PolytopeA = convexHull transpose matrix flowPolytopeVertices(QuiverA);

    for j from i  to (length RT - 1) do {
        -- Check if j is already in any of the sublists of AllIsoCases
        alreadyProcessed = any(AllIsoCases, sublist -> member(j, sublist));

        if not alreadyProcessed then {
            wb = RT#j;
            fb = incInverse(wb, K) + {-7, -11, 18, 7, 11, -18};
            QuiverB = toricQuiver(Q2D, fb);
            PolytopeB = convexHull transpose matrix flowPolytopeVertices(QuiverB);

            if areIsomorphic(PolytopeA, PolytopeB) == true then {
                -- Append j to the IsoCases for i
                AllIsoCases#i = append(AllIsoCases#i, j);
                print("Two isomorphic cases ", i, j);
            }
        }
    }
}

-- Use AllIsoCases as needed
-- Print only non-empty IsoCases
scan(AllIsoCases, IsoCases -> if #IsoCases > 0 then print IsoCases) 

for IsoCases in AllIsoCases do {
	if #IsoCases > 0 then 
	{
	print IsoCases;
	f << IsoCases << endl	
	}
}

 f << close
 


end
	
for i from 0 to (length RT-1) do  
	{	
	IsoCases = {i};
	--print(IsoCases);
	wa = RT#i;
	fa = incInverse(wa,K) + {-7,-11,18,7,11,-18};
	QuiverA = toricQuiver(Q2D, fa);	
	PolytopeA = convexHull transpose matrix flowPolytopeVertices( QuiverA );
	for j from i+1 to (length RT-1) do
		{
		wb = RT#j;
		fb = incInverse(wb,K) + {-7,-11,18,7,11,-18};
		QuiverB = toricQuiver(Q2D, fb);		
		PolytopeB = convexHull transpose matrix flowPolytopeVertices( QuiverB );
		if areIsomorphic(PolytopeA, PolytopeB ) == true then
			{		
			IsoCases = append(IsoCases, j);
			print("Two isomorphic cases ",i, j );
			}
		}			
	}	
	




AllData = {}
for i from 0 to (length RT-1) do  
	{	
		w = RT#i;
		fo = incInverse(w,K);
		f = fo + {-7,-11,18,7,11,-18};
		--f = {1,1,1,1,1,1};
		print("Case", i);
		print("Weight =", w);
		print("Flow = ",f);
				
		----
		Quiver2D = toricQuiver(Q2D, f);
		print( Quiver2D );	
		-- print (basisForFlowPolytope(Quiver2D)); -- it seems trouble some.	
		print( flowPolytopeVertices( Quiver2D ) );
		
		
		-- WARNING: basisForFlowPolytope(Quiver2D) seems to give the wrong answer when varying the flow.
		-- see how it does not change Polytope2D = convexHull transpose basisForFlowPolytope( Quiver2D );
		-- WARNING: flowPolytopeVertices seems to have problems when some of the weights are equal to zero.	
		
		-- Combinatorial information
		
		
		Polytope2D = convexHull transpose matrix flowPolytopeVertices( Quiver2D );		
		print("Ambient dimension polytope =", ambDim Polytope2D , "Compactness =", isCompact Polytope2D);		
		print( "Number of vertices =", nVertices Polytope2D );
		print("Number Lattice points =", length (latticePoints Polytope2D) );
		
		-- Geometric Information
		
		T = normalToricVariety transpose ( basisForFlowPolytope( Quiver2D) );
		dimT = dim T;		
		print("dimension Toric Varieties =", dimT); 
		print("Projective =", isProjective T); 
		print("Smooth =", isSmooth T); 
		print("Fano =", isFano T); 
		print("Projective space ", areIsomorphic(simplex dimT,Polytope2D ));
		
		--Data = {i, w, f, nVertices Polytope2D, length (latticePoints Polytope2D), dimT, isFano T, areIsomorphic(simplex dimT,Polytope2D ) };
		
		Data = {i, w, nVertices Polytope2D, length (latticePoints Polytope2D), isFano T, areIsomorphic(simplex dimT,Polytope2D ), dim picardGroup T };
		
		AllData = append(AllData, Data)			
			
	}
	

print(AllData);	
	

end


-- the following function calculates a random integer linear combination
-- of the columns of a matrix. 
RLC = N  -> 
    (
    n = numColumns(N); -- Get the length of the list.
    result = (random(1,5))*N_0; -- Initialize the result.
    
    for i from 1 to (n-1) do 
    	{
    	c = (random(1,5))*n; -- Generate random coefficients 
        result = result + c*N_i; 
        };
    result -- Return the result.
    )



M = quiverIncidenceMatrix(K) 

-- The columns are the primitive weights, and the number of rows -1 is the dimension
-- where the cone C(Q) exist. So this is the number of rays required to span a maximal
-- fan.

S = subsets((numColumns M),(numRows M)-1); 

dQ = (numColumns M) - (numRows M) +1 

print(length S)

ListV = {}; -- The number of integers
ListFano = {}; -- If the associated toric varieties is Fano

print(dQ)

for i from 0 to max(length(S)-1,2) do 
	{
		print("Case", i);	
		-- select a random element of the camber
		w = entries RLC(M_(S#i));
		
		--print(M_(S#i));
		--f = incInverse(w,K);
		--print(i,w,f); 
		
		-- construct the polytope
		MT = transpose matrix flowPolytope(w,K);
		CT = convexHull MT;
		if dim CT == dQ then -- testing if polytope full dimensional
			(		
			-- construct the toric variety associated to that chamber
			print("Weight =", w);
			T = normalToricVariety MT;
			nv = nVertices (convexHull MT);
			print("# vertices", nv );
			ListV =  append(ListV,nv);
			
			-- next properties toric varieties associated to each chamber			
			print("Smooth =", isSmooth T);
			print("Projective =", isProjective T);
			print("Fano =", isFano T);
			ListFano = append(ListFano, isFano T);
			print("Dim =", dim T);
			
			PT = convexHull transpose matrix flowPolytope(w,K);
			print(transpose matrix flowPolytope(w,K));
			print("Projective space ", areIsomorphic(simplex dim T,PT));
			print(" ");
			)
	}

print(ListV)

print(ListFano)

stop()





