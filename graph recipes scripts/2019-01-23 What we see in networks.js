// Init
var settings = {}
settings.save = false

// Compute Louvain modularity (Graphology)
louvain.assign(g)

// Count Louvain classes
var louvainClasses = {}
g.nodes().forEach(function(nid){
	var louvainClass = g.getNodeAttribute(nid, "community")
	louvainClasses[louvainClass] = true
})

// Rename Louvain classes for clarity
var alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
Object.keys(louvainClasses).forEach(function(louvainClass, i){
	louvainClasses[louvainClass] = i
})
g.nodes().forEach(function(nid){
	g.setNodeAttribute(nid, "community", alphabet[louvainClasses[g.getNodeAttribute(nid, "community")]])
})


// Run a k-means initiated with random values
var kmeans_iter = 100
var kmeans_k = Object.keys(louvainClasses).length
var means = []
for (let i=0; i<kmeans_k; i++) {
	var ni = Math.floor(Math.random() * g.order)
	while (means.some(function(m){m==ni})) {
		ni = Math.floor(Math.random() * g.order)
	}
	means.push(ni)
}
means = means.map(function(ni){ return g.getNodeAttributes(g.nodes()[ni]) })
means = kmeans(means, 100)

// Attribute kmeans class to the closest centroid
g.nodes().forEach(function(nid){
	var n = g.getNodeAttributes(nid)
	var closest = 0
	var d2min = dist2(n, means[closest]) 
	means.forEach(function(m, i){
		if (i>0) {
			var d2 = dist2(n, m)
			if (d2 < d2min) {
				d2min = d2
				closest = i
			}
		}
	})
	n.kclass = closest
})
g.nodes().forEach(function(nid){
	g.setNodeAttribute(nid, "kclass", alphabet[g.getNodeAttribute(nid, "kclass")])
})

// Random class
g.nodes().forEach(function(nid){
	g.setNodeAttribute(nid, "randclass", alphabet[Math.floor(Math.random()*kmeans_k)])
})

// Compute the Jaccard indexes of pairs
var jaccardPairs = [
	{att1:"community", att2:"kclass"},
	{att1:"randclass", att2:"kclass"},
	{att1:"randclass", att2:"community"}
]
jaccardPairs.forEach(function(d){
	d.jaccardIndex = computeJaccard(d.att1, d.att2)
})

// Draw the networks
initNetworkDrawing()
d3.select('#playground').append('h4').text('Jaccard indexes')
d3.select('#playground').append('p').append('em').text('This number between 0 and 1 measures at which point the two partitions are similar. It is defined as the Jaccard index of the pairs of nodes that are in the same cluster.')
jaccardPairs.forEach(function(d){
	d3.select('#playground').append('p').text('Index of '+d.att1+'-'+d.att2+': '+d.jaccardIndex)
})
d3.select('#playground').append('h4').text('Network colored by Louvain modularity clustering')
drawNetwork({nodecolor: 'attributeQuali', attribute: 'community', label: 'community', fileSuffix: ' Louvain Class'})
d3.select('#playground').append('h4').text('Network colored by k-means (k = number of Louvain clusters)')
drawNetwork({nodecolor: 'attributeQuali', attribute: 'kclass', label: 'kclass', fileSuffix: ' K Class'})
d3.select('#playground').append('h4').text('Network colored by random classes (same number of Louvain clusters)')
drawNetwork({nodecolor: 'attributeQuali', attribute: 'randclass', label: 'randclass', fileSuffix: ' Random Class'})





/// FUNCTIONS
function computeJaccard(att1, att2) {
	var union = 0
	var intersection = 0
	g.nodes().forEach(function(n1id){
		g.nodes().forEach(function(n2id){
			var att1Same = g.getNodeAttribute(n1id, att1) == g.getNodeAttribute(n2id, att1)
			var att2Same = g.getNodeAttribute(n1id, att2) == g.getNodeAttribute(n2id, att2)
			if (att1Same || att2Same) {
				union++
				if (att1Same && att2Same) {
					intersection++
				}
			}
		})
	})
	return intersection / union
}

function initMatrix(size, defaultValue) {
	var i
	var j
	var result = []
	for (i=0; i<size; i++) {
		var row = []
		for (j=0; j<size; j++) {
			row.push(defaultValue)
		}
		result.push(row)
	}
	return result
}

function rayleighQ(M, x) {
	return numeric.dot(numeric.dot(x, M), x) / numeric.dot(x, x)
}

function jacobi(a, it_max) {

	// Adapted from https://people.sc.fsu.edu/~jburkardt/py_src/jacobi_eigenvalue/jacobi_eigenvalue.py
  
  var c, d, g, h, i, j, k, l, m, n, p, q, s, t, v, w, bw, zw
	var gapq, it_num, rot_num, tau, term, termp, termq, theta, thresh
	
	n = a.length // matrix size
  v = initMatrix(n, 0)
  d = a.map(function(){return 0})
  
  for (j=0; j<n; j++) {
  	for (i=0; i<n; i++) {
  		v[i][j] = 0
  	}
		v[j][j] = 1
  }

	for (i=0; i<n; i++) {
		d[i] = a[i][i]
	}

	bw = d.slice(0)
	zw = a.map(function(){return 0})
	w = a.map(function(){return 0})

	for (i=0; i<n; i++) {
		bw[i] = d[i]
	}

  it_num = 0
  rot_num = 0

  while ( it_num < it_max ) {

    it_num++

		if (it_num%10000 == 0)
			console.log('Jacobi iteration ' + it_num + ' | state ' + state)
		
		//  The convergence threshold is based on the size of the elements in
		//  the strict upper triangle of the matrix.
    
    thresh = 0
    for (j=0; j<n; j++) {
	  	for (i=0; i<j; i++) {
	  		thresh += Math.pow(a[i][j], 2)
	  	}
	  }

	  thresh = Math.sqrt(thresh) / (4 * n)

	  if (thresh == 0) {
	  	break
	  }

	  for (p=0; p<n; p++) {
	  	for (q=p+1; q<n; q++) {
	  		gapq = 10 * Math.abs( a[p][q] )
        termp = gapq + Math.abs( d[p] )
        termq = gapq + Math.abs( d[q] )

				if ( 4 < it_num && termp == Math.abs( d[p] ) && termq == Math.abs( d[q] ) ) {
					//  Annihilate tiny offdiagonal elements.
					a[p][q] = 0

				} else if(thresh <= Math.abs( a[p][q] )) {
					//  Otherwise, apply a rotation.
					h = d[q] - d[p]
          term = Math.abs( h ) + gapq

          if (term == Math.abs(h)) {
          	t = a[p][q] / h
          } else {
          	theta = 0.5 * h / a[p][q]
						t = 1 / ( Math.abs( theta ) + Math.sqrt( 1 + theta * theta ) )
            if ( theta < 0 ){
              t = - t
            }
          }
          c = 1 / Math.sqrt( 1 + t * t )
          s = t * c
          tau = s / ( 1 + c )
          h = t * a[p][q]

					//  Accumulate corrections to diagonal elements.

					zw[p] = zw[p] - h                  
          zw[q] = zw[q] + h
          d[p] = d[p] - h
          d[q] = d[q] + h

          a[p][q] = 0

          //  Rotate, using information from the upper triangle of A only.

          for (j=0; j<p; j++) {
          	g = a[j][p]
          	h = a[j][q]
          	a[j][p] = g - s * ( h + g * tau )
            a[j][q] = h + s * ( g - h * tau )
          }

          for (j=p+1; j<q; j++) {
          	g = a[p][j]
          	h = a[j][q]
          	a[p][j] = g - s * ( h + g * tau )
            a[j][q] = h + s * ( g - h * tau )
          }

          for (j=q+1; j<n; j++) {
          	g = a[p][j]
          	h = a[q][j]
          	a[p][j] = g - s * ( h + g * tau )
            a[q][j] = h + s * ( g - h * tau )
          }

					//  Accumulate information in the eigenvector matrix.

					for (j=0; j<n; j++) {
						g = v[j][p]
            h = v[j][q]
            v[j][p] = g - s * ( h + g * tau )
            v[j][q] = h + s * ( g - h * tau )
					}

					rot_num++
				}
	  	}
	  }

	  for (i=0; i<n; i++) {
      bw[i] = bw[i] + zw[i]
      d[i] = bw[i]
      zw[i] = 0	  	
	  }
	}
	
	//  Restore upper triangle of input matrix.
  
	for (j=0; j<n; j++) {
  	for (i=0; i<j; i++) {
  		a[i][j] = a[j][i]
  	}
  }

	//  Ascending sort the eigenvalues and eigenvectors.

	for (k=0; k<n-1; k++) {
		m = k
		for (l=k+1; l<n; l++) {
			if ( d[l] < d[m] ) {
        m = l
			}
		}
		if ( l != m ) {

			t    = d[m]
      d[m] = d[k]
      d[k] = t

      for (i=0; i<n; i++) {
      	w[i]   = v[i][m]
        v[i][m] = v[i][k]
        v[i][k] = w[i]
      }

		}
	}

	// console.log('Jacobi data: v', v, 'd', d, 'it_num', it_num, 'rot_num', rot_num)

	// Format result
	
	var eigenvectors = {}
	d.forEach(function(eigval, i){
		var vectors = eigenvectors[eigval] || []
		var eigenvector = v.map(function(row){
			return row[i]
		})
		vectors.push(eigenvector)
		eigenvectors[eigval] = vectors
	})
	var eigenvalues = d.slice(0).sort(function(a, b){
		return a - b
	})
	return {eigenvalues:eigenvalues, eigenvectors:eigenvectors}

	return {eigenvalues:eigenvalues, eigenvectors:eigenvectors}
}


function initNetworkDrawing() {
	// Canvas size
	settings.width =  1000
	settings.height = 1000
	settings.offset = 20 // Margin

	// Drawing nodes, labels and edges
	settings.display_label = true
	settings.node_size = 4
	settings.font_size = 14
	settings.font_family = 'Open Sans Condensed, sans-serif'
	settings.font_weight = 300
	settings.edge_color = 'rgba(200, 200, 200, 0.3)'

	// --- (end of settings)

	// Fix missing coordinates and/or colors
	addMissingVisualizationData()

	// Change the coordinates of the network to fit the canvas space
	rescaleGraphToGraphicSpace()

	// Init graphic space
	document.querySelector('#playground').innerHTML = ''
}

function rescaleGraphToGraphicSpace() {

  // Flip Y because canvas drawing differs from Maths
  g.nodes().forEach(function(nid){
  	var n = g.getNodeAttributes(nid)
    n.y = -n.y
  })

  // General barycenter resize
  var xbarycenter = 0
  var ybarycenter = 0
  var wtotal = 0
  var dx
  var dy
  var ratio

  g.nodes().forEach(function(nid){
  	var n = g.getNodeAttributes(nid)
    // We use node size as weight (default to 1)
    n.size = n.size || 1
    xbarycenter += n.size * n.x
    ybarycenter += n.size * n.y
    wtotal += n.size
  })
  xbarycenter /= wtotal
  ybarycenter /= wtotal

  var dmax = 0 // Maximal distance from barycenter
  g.nodes().forEach(function(nid){
  	var n = g.getNodeAttributes(nid)
    var d = Math.sqrt( Math.pow(n.x - xbarycenter, 2) + Math.pow(n.y - ybarycenter, 2) )
    dmax = Math.max(dmax, d)
  })

  ratio = ( Math.min(settings.width, settings.height) - 2 * settings.offset ) / (2 * dmax)

  // Initial resize
  g.nodes().forEach(function(nid){
  	var n = g.getNodeAttributes(nid)
    n.x = settings.width / 2 + (n.x - xbarycenter) * ratio
    n.y = settings.height / 2 + (n.y - ybarycenter) * ratio
    n.size *= ratio
  })
}

function addMissingVisualizationData() {
  var colorIssues = 0
  var coordinateIssues = 0
  g.nodes().forEach(function(nid){
    var n = g.getNodeAttributes(nid)
    if (!isNumeric(n.x) || !isNumeric(n.y)) {
      var c = getRandomCoordinates()
      n.x = c[0]
      n.y = c[1]
      coordinateIssues++
    }
    if (!isNumeric(n.size)) {
      n.size = 1
    }
    if (n.color == undefined) {
      n.color = '#665'
      colorIssues++
    }
  })

  if (coordinateIssues > 0) {
    alert('Note: '+coordinateIssues+' nodes had coordinate issues. We carelessly fixed them.')
  }

  function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n)
  }
  
  function getRandomCoordinates() {
    var candidates
    var d2 = Infinity
    while (d2 > 1) {
      candidates = [2 * Math.random() - 1, 2 * Math.random() - 1]
      d2 = candidates[0] * candidates[0] + candidates[1] * candidates[1]
    }
    var heuristicRatio = 5 * Math.sqrt(g.order)
    return candidates.map(function(d){return d * heuristicRatio})
  }
}

function drawNetwork(opt) {
	opt = opt || {}

	// Create the canvas
	var div = document.createElement("div")
	div.setAttribute('style', 'width: ' + settings.width + 'px; height: ' + settings.height + 'px;')
	document.querySelector('#playground').appendChild(div)

	var canvas = document.createElement("canvas")
	canvas.setAttribute('width', settings.width)
	canvas.setAttribute('height', settings.height)
	div.appendChild(canvas)

	var ctx = canvas.getContext("2d")

	// Paint a white background
	ctx.beginPath()
	ctx.rect(0, 0, settings.width, settings.height)
	ctx.fillStyle="white"
	ctx.fill()
	ctx.closePath()

	// Draw each edge
	g.edges().forEach(function(eid){
		var ns = g.getNodeAttributes(g.source(eid))
		var nt = g.getNodeAttributes(g.target(eid))
		var e = g.getEdgeAttributes(eid)

		var color
		if (opt.edgecolor == 'original') {
			color = e.color
		} else {
			color = settings.edge_color
		}

	  ctx.beginPath()
	  ctx.lineCap="round"
	  ctx.lineJoin="round"
	  ctx.strokeStyle = color
	  ctx.fillStyle = 'rgba(0, 0, 0, 0)';
	  ctx.lineWidth = settings.edge_thickness
	  ctx.moveTo(ns.x, ns.y)
	  ctx.lineTo(nt.x, nt.y)
	  ctx.stroke()
	  ctx.closePath()
	})

	// Reset nodes color according to options
	if (opt.nodecolor == 'original') {
		// Nothing to do
	} else if (opt.nodecolor == 'attribute') {
		var ext = d3.extent(g.nodes().map(function(nid){return g.getNodeAttribute(nid, opt.attribute)}))
		var nodeColor = d3.scaleLinear()
			.domain([ext[0], ext[1]])
			.range(['orange', 'purple'])
		g.nodes().forEach(function(nid, i){
			var n = g.getNodeAttributes(nid)
			n.color = nodeColor(n[opt.attribute])
		})
	} else if (opt.nodecolor == 'attributeQuali') {
		var dcolor = g.nodes().map(function(nid){return g.getNodeAttribute(nid, opt.attribute)})
		var nodeColor = d3.scaleOrdinal(d3.schemeCategory10)
			.domain(dcolor)
			
		g.nodes().forEach(function(nid, i){
			var n = g.getNodeAttributes(nid)
			n.color = nodeColor(n[opt.attribute])
		})
	} else if (opt.nodecolor == 'attributeOppositive') {
		var ext = d3.extent(g.nodes().map(function(nid){return g.getNodeAttribute(nid, opt.attribute)}))
		var nodeColor = d3.scaleLinear()
			.domain([ext[0], 0, ext[1]])
			.range(['red', '#EEEEEE', 'green'])
		g.nodes().forEach(function(nid, i){
			var n = g.getNodeAttributes(nid)
			n.color = nodeColor(n[opt.attribute])
		})
	} else {
		g.nodes().forEach(function(nid, i){
			var n = g.getNodeAttributes(nid)
			n.color = "#666"
		})
	}

	// Draw each node
	g.nodes().forEach(function(nid){
		var n = g.getNodeAttributes(nid)

	  ctx.lineCap="round"
	  ctx.lineJoin="round"

	  if (settings.display_label) {
	    ctx.font = settings.font_weight + " " + settings.font_size+"px "+settings.font_family;
	    ctx.lineWidth = 4
	    ctx.fillStyle = '#FFFFFF'
	    ctx.strokeStyle = '#FFFFFF'
	    var label = n[opt.label]
	    if (opt.round) {
	    	label = (+label).toPrecision(3)
	    }
	    ctx.fillText(
	      label
	    , n.x + settings.node_size * 1.4
	    , n.y + 0.3 * settings.font_size
	    )
	    ctx.strokeText(
	      label
	    , n.x + settings.node_size * 1.4
	    , n.y + 0.3 * settings.font_size
	    )
	    ctx.lineWidth = 0
	    ctx.fillStyle = n.color
	    ctx.fillText(
	      label
	    , n.x + settings.node_size * 1.4
	    , n.y + 0.3 * settings.font_size
	    )
	  }

	  ctx.beginPath()
	  ctx.arc(n.x, n.y, settings.node_size, 0, 2 * Math.PI, false)
	  ctx.lineWidth = 0
	  ctx.fillStyle = n.color
	  ctx.shadowColor = 'transparent'
	  ctx.fill()
	})

	// Save if needed
	if (settings.save) {
		saveCanvas(canvas, store.get('graphname') + (opt.fileSuffix || '') )
	}
}

function drawScatterplot(data, opt) {
	opt = opt || {}

	var container = d3.select('#playground').append('div')

	var margin = {top: 24, right: 24, bottom: 32, left: 32}
	var width = opt.width - margin.left - margin.right
	var height = opt.height  - margin.top - margin.bottom

	var x = d3.scaleLinear().range([0, width])
	var y = d3.scaleLinear().range([height, 0])

	var xAxis = d3.axisBottom()
	    .scale(x)

	var yAxis = d3.axisLeft()
	    .scale(y)

	var svg = container.append("svg")
	    .attr("width", width + margin.left + margin.right)
	    .attr("height", height + margin.top + margin.bottom)
	  .append("g")
	    .attr("transform", 
	        "translate(" + margin.left + "," + margin.top + ")")

	x.domain(d3.extent(data.map(function(d){return d[opt.attributeX]})))
	y.domain(d3.extent(data.map(function(d){return d[opt.attributeY]})))

	svg.selectAll("dot")
	      .data(data)
	    .enter().append("circle")
	      .attr("r", 3)
	      .attr("cx", function(d) { return x(d[opt.attributeX]); })
	      .attr("cy", function(d) { return y(d[opt.attributeY]); })
	      .attr('fill', 'rgba(0, 0, 0, 0.05)')

	// Add the X Axis
	svg.append("g")
	    .attr("transform", "translate(0," + height + ")")
	    .call(d3.axisBottom(x))
	  .append("text")
      .attr("x", width)
      .attr("y", -6)
      .style("text-anchor", "end")
      .attr('fill', '#000')
      .attr('font-family', 'Roboto,Helvetica Neue,sans-serif')
      .attr('font-size', '12px')
      .text(opt.labelX || '')

	// Add the Y Axis
	svg.append("g")
	    .call(d3.axisLeft(y))
	  .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .attr('fill', '#000')
      .attr('font-family', 'Roboto,Helvetica Neue,sans-serif')
      .attr('font-size', '12px')
      .text(opt.labelY || '')

	if (settings.save) {
		saveSVG(svg, store.get('graphname') + ' Scatterplot')
	}
}

function correlation(a1, a2) {
	// Covariance
	var cov = d3.mean(a1.map(function(a1_i,i){return a1_i*a2[i]})) - d3.mean(a1)*d3.mean(a2)
	return cov / (d3.deviation(a1) * d3.deviation(a2))
}

function saveSVG(svg, name) {
	// Download SVG
	var svgFileContent = []
	svgFileContent.push('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="'+settings.width+'" height="'+settings.height+'" viewBox="0 0 '+settings.width+' '+settings.height+'">')
	svgFileContent.push(svg.html())
	svgFileContent.push('</svg>')

	var blob = new Blob(svgFileContent, {type: "image/svg+xml;charset=utf-8"})
	saveAs(blob, name + ".svg")
}

function saveCanvas(canvas, name) {
  canvas.toBlob(function(blob) {
    saveAs(blob, name + ".png");
  });
}

function kmeans(means, iter) {
	while (iter-->0) {
		var barycenters = means.map(function(){ return {xsum:0, ysum:0, total:0} })
		// Find distribution
		g.nodes().forEach(function(nid, i){
			var n = g.getNodeAttributes(nid)
			var minD = Infinity
			var cj
			means.map(function(c, j){
				var d = Math.sqrt(Math.pow(c.x - n.x, 2) + Math.pow(c.y - n.y, 2))
				if (d < minD) {
					minD = d
					cj = j
				}
			})
			barycenters[cj].xsum += n.x
			barycenters[cj].ysum += n.y
			barycenters[cj].total++
		})
		// Reattribute centroids
		means = barycenters.map(function(b){ return {x: b.xsum/b.total, y: b.ysum/b.total} })
	}
	return means
}

function dist2_id(n1id, n2id) {
	var x1 = g.getNodeAttribute(n1id, "x")
	var y1 = g.getNodeAttribute(n1id, "y")
	var x2 = g.getNodeAttribute(n2id, "x")
	var y2 = g.getNodeAttribute(n2id, "y")
	return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
}

function dist2(n1, n2) {
	return (n1.x-n2.x)*(n1.x-n2.x) + (n1.y-n2.y)*(n1.y-n2.y)
}

function dist_id(n1id, n2id) { return Math.sqrt(dist2_id(n1id, n2id)) }

function dist(n1, n2) { return Math.sqrt(dist2(n1, n2)) }