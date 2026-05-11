function getEdgeMutationCount(d) {
  if (!d || d.changes === undefined || d.changes === null) return 0;

  var changes = String(d.changes).replace(/^\s+|\s+$/g, '');
  var count;

  if (changes === '') return 0;
  if (/^\d+(\.\d+)?$/.test(changes)) {
    count = Math.round(Number(changes));
    if (count > 0) return count;
    return 0;
  }

  return 1;
}

function getEdgeMutationSegments() {
  var segments = [];
  var count;

  for (var i = 0; i < linkList.length; i++) {
    count = getEdgeMutationCount(linkList[i]);
    for (var j = 0; j < count; j++) {
      segments.push({ link: linkList[i], index: j, count: count });
    }
  }

  return segments;
}

function updateEdgeMutationElements() {
  var labelOffset = 10;
  var segmentSpacing = 7;
  var segmentLength = 8;

  svg.selectAll('.edge-mutation-label')
     .attr('x', function(d) {
       var dx = d.target.x - d.source.x;
       var dy = d.target.y - d.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       return (d.source.x + d.target.x) / 2 - (dy / len) * labelOffset;
     })
     .attr('y', function(d) {
       var dx = d.target.x - d.source.x;
       var dy = d.target.y - d.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       return (d.source.y + d.target.y) / 2 + (dx / len) * labelOffset;
     });

  svg.selectAll('.edge-mutation-segment')
     .attr('x1', function(d) {
       var dx = d.link.target.x - d.link.source.x;
       var dy = d.link.target.y - d.link.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       var ux = dx / len;
       var uy = dy / len;
       var nx = -uy;
       var shift = (d.index - (d.count - 1) / 2) * segmentSpacing;
       var cx = (d.link.source.x + d.link.target.x) / 2 + ux * shift;
       return cx - nx * segmentLength / 2;
     })
     .attr('y1', function(d) {
       var dx = d.link.target.x - d.link.source.x;
       var dy = d.link.target.y - d.link.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       var ux = dx / len;
       var uy = dy / len;
       var ny = ux;
       var shift = (d.index - (d.count - 1) / 2) * segmentSpacing;
       var cy = (d.link.source.y + d.link.target.y) / 2 + uy * shift;
       return cy - ny * segmentLength / 2;
     })
     .attr('x2', function(d) {
       var dx = d.link.target.x - d.link.source.x;
       var dy = d.link.target.y - d.link.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       var ux = dx / len;
       var uy = dy / len;
       var nx = -uy;
       var shift = (d.index - (d.count - 1) / 2) * segmentSpacing;
       var cx = (d.link.source.x + d.link.target.x) / 2 + ux * shift;
       return cx + nx * segmentLength / 2;
     })
     .attr('y2', function(d) {
       var dx = d.link.target.x - d.link.source.x;
       var dy = d.link.target.y - d.link.source.y;
       var len = Math.sqrt(dx * dx + dy * dy) || 1;
       var ux = dx / len;
       var uy = dy / len;
       var ny = ux;
       var shift = (d.index - (d.count - 1) / 2) * segmentSpacing;
       var cy = (d.link.source.y + d.link.target.y) / 2 + uy * shift;
       return cy + ny * segmentLength / 2;
     });
}

function updateSVG() {
  link = svg.selectAll('.link').remove();
  svg.selectAll('.edge-mutation-label').remove();
  svg.selectAll('.edge-mutation-segment').remove();

  link = svg.selectAll('.link').data(linkList);
  link.enter().append('line')
              .attr('class', 'link')
              .attr('id', function(d) { return d.name; } )
              .on('mouseover', function(){ d3.select(this).style({'stroke': '#FF011B', 'stroke-width': linewidth*3}); })
              .on('mouseout', function(){ d3.select(this).style({'stroke': '#000000', 'stroke-width': linewidth});})
              .on('click', clickLink)
              .append('title')
              .text(function(d, i) { return d.changes; });
              
  link.style('stroke-width', linewidth ).style('stroke', '#000000');           
  link.exit().remove();

  if (edgeMutationDisplayMode === 'number') {
    svg.selectAll('.edge-mutation-label')
       .data(linkList.filter(function(d) { return getEdgeMutationCount(d) > 0; }))
       .enter()
       .append('text')
       .attr('class', 'edge-mutation-label')
       .style('font-size', '11px')
       .style('font-family', 'Helvetica, Arial, sans-serif')
       .style('font-weight', 'bold')
       .style('text-anchor', 'middle')
       .style('dominant-baseline', 'central')
       .style('paint-order', 'stroke')
       .style('stroke', '#ffffff')
       .style('stroke-width', '3px')
       .style('fill', '#000000')
       .style('pointer-events', 'none')
       .text(function(d) { return '(' + getEdgeMutationCount(d) + ')'; });
  } else if (edgeMutationDisplayMode === 'segments') {
    svg.selectAll('.edge-mutation-segment')
       .data(getEdgeMutationSegments())
       .enter()
       .append('line')
       .attr('class', 'edge-mutation-segment')
       .style('stroke', '#000000')
       .style('stroke-width', Math.max(1, Number(linewidth)))
       .style('stroke-linecap', 'round')
       .style('pointer-events', 'none');
  }

  node = svg.selectAll('.node').remove();
  node = svg.selectAll('.node').data(nodeList);
  node.enter().append('g')
              .attr('class', 'node')
              .attr('id', function(d) { return d.name; } )
              .on('click', clickNode)
              .call(drag)
              .insert('circle')
              .attr('r', function(d){ return d.radius;})
              .on('mouseover', function(){ d3.select(this).style({'stroke': '#FF011B', 'stroke-width': linewidth*4 }); })
              .on('mouseout', function(){d3.select(this).style({'stroke': '#000000', 'stroke-width': linewidth*2 });})
              .append('title')
              .text(function(d, i) { return d.name; });
  node.style('stroke-width', function(){ return linewidth*2; }).style('stroke', '#000000').style('fill', 'white');
  node.exit().remove();
              
  path = node.selectAll('path').remove();
  path = node.selectAll('path').data(function(d, i) { return pie(d.proportions); });
  path.enter()
      /*
       * Using the construct below, the arcs are inserted before the circle of the
       * node. It allows for some interactions (e.g., coloring the whole circle, etc,
       * and avoids the css use of 'path { pointer-events: none; }' to allow the
       * node (which is below) to capture 'mouseover' and 'mouse out' events
       */
      .append('path')
      .attr('d', sector)
      .style('fill', function(d) { 
        if (d.data.pattern === 'none' ) return d.data.color;
        else return d.data.pattern;
      });
      
      /*
       * These two commands may be used to implement
       * a stroke between arcs in the pie
       * .style('stroke-width', '0')
       * .style('stroke', 'none');
       */
       if (outlinenodes) {
         path.style('stroke-width', linewidth/2 ).style('stroke', '#000000');
       } else {
         path.style('stroke-width', '0').style('stroke', 'none');
       }
  path.exit().remove();
  updateEdgeMutationElements();
  
  force.nodes(nodeList).links(linkList).start();
}
