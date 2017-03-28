/**
 * Created by asherkhb on 3/7/17.
 */

var graph;
var data;

var limit = 1;

var words = [];
var linked = [];

var group2color = {'1': 'coral', '2': 'aquamarine', '3': 'lavender'};


function make_graph(el) {

    // Add and remove elements on the graph object
    this.addNode = function (id, r, c) {
        nodes.push({"id":id, "r":r, "color":c});
        update();
    };

    this.removeNode = function (id) {
        var i = 0;
        var n = findNode(id);
        while (i < links.length) {
            if ((links[i]["source"] == n)||(links[i]["target"] == n))
            {
                links.splice(i,1);
            }
            else i++;
        }
        nodes.splice(findNodeIndex(id),1);
        update();
    };

    this.removeLink = function (source,target){
        for(var i=0;i<links.length;i++)
        {
            if(links[i].source.id == source && links[i].target.id == target)
            {
                links.splice(i,1);
                break;
            }
        }
        update();
    };

    this.removeallLinks = function(){
        links.splice(0,links.length);
        update();
    };

    this.removeAllNodes = function(){
        nodes.splice(0,links.length);
        update();
    };

    this.addLink = function (source, target, value) {
        links.push({"source":findNode(source),"target":findNode(target),"value":value});
        update();
    };

    var findNode = function(id) {
        for (var i in nodes) {
            if (nodes[i]["id"] === id) return nodes[i];};
    };

    var findNodeIndex = function(id) {
        for (var i=0;i<nodes.length;i++) {
            if (nodes[i].id==id){
                return i;
            }
        };
    };

    // set up the D3 visualisation in the specified element
    var w = 600,
        h = 600;
    var vis = d3.select(el)
        .append("svg:svg")
        .attr("width", w)
        .attr("height", h)
        .attr("id","svg")
        .attr("pointer-events", "all")
        //.attr("viewBox","0 0 "+w+" "+h)
        .attr("perserveAspectRatio","xMinYMid")
        .append("svg:g");

    var force = d3.layout.force();

    var nodes = force.nodes(),
        links = force.links();

    var update = function () {
        var link = vis.selectAll("line")
            .data(links, function(d) {
                return d.source.id + "-" + d.target.id;
                });

        link.enter().append("line")
            .attr("id",function(d){return d.source.id + "-" + d.target.id;})
            .attr("opacity", 1)
            .attr("stroke-width", 1)
            .attr("stroke", "blue")
            .transition().duration(500)
            .attr("opacity", 0.3)
            .attr("stroke-width", function(d){return (Math.abs(d.value)+1)*1.5})
            .attr("stroke", function(d){
                if (d.value>0){return "teal"}else{return "pink"}})
            .attr("class","link");
        link.append("title")
        .text(function(d){return d.value;});

        link.exit().transition().delay(500).remove();





        var node = vis.selectAll("g.node")
            .data(nodes, function(d) {
                return d.id;}); //.on('mouseover', showAaInfo);

        var nodeEnter = node.enter().append("g")
            .attr("class", "node")
            .call(force.drag);

        nodeEnter.append("svg:circle")
        .attr("r", function(d){return d.r*3})
        .attr("id",function(d) { return "Node_"+d.id;})
        .attr("fill", function(d){return d.color})
        .attr("class","nodeStrokeClass");

        nodeEnter.append("svg:text")
        .attr("text-anchor", "middle")
        .attr("class","textClass")
        .attr("id",function(d) { return "Node_"+d.id;})
        .attr("fill", "black")
        .attr("font-size", function(d){return d.r*6})
        .attr("text-align", "center")
        //.attr("dx", function(d){return -d.r*1.2})
        .attr("dy", function(d){return d.r*1.9})
        .text( function(d){return d.id;});

        node.exit().remove();
        force.friction(0.9);
        force.on("tick", function() {



            node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y         + ")"; });

            link.attr("x1", function(d) { return d.source.x; })
              .attr("y1", function(d) { return d.source.y; })
              .attr("x2", function(d) { return d.target.x; })
              .attr("y2", function(d) { return d.target.y; });


        });

        // Restart the force layout.
        force
        .gravity(.2)
        .friction(0.7)
        .charge(-1200)
        .distance(100)
        .linkDistance(100)
        .size([w, h])
        .start();
    };


    // Make it all go
    update();

}

function initialize() {
    // Build graph.
    graph = new make_graph("#graph");
    data = [];
	d3.csv("embedtest.csv", function(curr){
        // Save data from CSV.
        data = curr;

        // Add Nodes.
        for(var i=0; i<data.length; i++) {
            var word = data[i][""];
            words[i] = word;

            var radius = data[i][word] / 1000;  // Radius = count / 1000. TODO: Intelligent radius calculation.
            var color = group2color['1'];  // Color from list. TODO: Intelligent color determination.
            graph.addNode(word, radius, color);
        }
        var x = words.length;
        for (var j=0; j<x; j++) {
            for (var k=j; k<x; k++) {
                var word_a = words[j];
                var word_b = words[k];

                var similarity = data[k][word_a];
                if (similarity >= limit) {
                    graph.addLink(word_a, word_b, similarity);
                    linked.push(word_a + "||" + word_b); // Take note that values are currently linked.
                }
            }
        }
    });

    // Set up slider.
    var sliderText = [
        {"label":"SIMILARITY SCORE", "x":75, "y":-3, "fs":"1em", "fc":"gainsboro"},
        {"label":"0", "x":0, "y":20, "fs":"0.9em", "fc":"black"},
        {"label":".5", "x":140, "y":20, "fs":"0.9em", "fc":"black"},
        {"label":"1", "x":285, "y":20, "fs":"0.9em", "fc":"black"}
    ];
    var sl = d3.select("#sliderInfo").append("svg")
                        .attr("height", 35)
                        .attr("width", 300)
                        .attr("color", "ligthgray");

    sl.selectAll("text").data(sliderText).enter().append("svg:text")
                        .attr("class", "sliderText")
                        .text(function(d){return d.label})
                        .attr("x", function(d){return d.x} )
                        .attr("y", function(d){return d.y+15})
                        .attr("font-size", function(d){return d.fs})
                        .attr("fill", function(d){return d.fc});
}

function updateLinks(new_limit) {
    var x = words.length;
    for (var i=0; i<x; i++) {
        for (var j=i; j<x; j++) {
            var word_a = words[i];
            var word_b = words[j];
            var similarity = data[j][word_a];
            var current_link = linked.indexOf(word_a + "||" + word_b);

            if (similarity >= new_limit) {
                if (current_link < 0) {
                    graph.addLink(word_a, word_b, similarity);
                    linked.push(word_a + "||" + word_b);
                }
            } else {
                if (current_link < 0) {
                    graph.removeLink(word_a, word_b);
                    linked.splice(current_link, 1);
                }
            }
        }
    }
}

function changeLimit(){
	limit = d3.select("#limit").property("value") * 0.01;
	updateLinks(parseFloat(limit));
}