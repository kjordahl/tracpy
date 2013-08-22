window.onload = function () {

    var cloudmadeUrl = 'http://{s}.tile.cloudmade.com/4fd891040a8a4ecb805c388019e23d46/64082/256/{z}/{x}/{y}.png';

    var basemap = new L.TileLayer(cloudmadeUrl, {maxZoom: 18});
    var latlng = new L.LatLng(27, -94);
    var delay = 40;		// animation delay (larger is slower)
    var Npoints = 151;		// number of points per track
    var isRunning = true;
    var plot_tracks = true;
    var plot_markers = true;

    function showTimeStep(j) {
	if (plot_tracks) {
	    tracks.eachLayer(function (layer) {
		showTrackStep(layer, j);
	    });
	}
	if (plot_markers) {
	    markers.eachLayer(function (layer) {
		showMarkerStep(layer, j);
	    });
	}
	if (isRunning) {
	    nextj = (j + 1) % Npoints;
	    setTimeout(function(i) {
		return function() {
		    showTimeStep(i);
		}
	    }(nextj), delay);
	}
    }

    function animateLines(data) {
	showTimeStep(0);
    }

    function showTrackStep(track, j) {
	pt = tracklines[track.track_id].geometry.coordinates[j];
	if (pt !== undefined) {
	    if (j === 0) {
		track.setLatLngs([[pt[1], pt[0]]]);
	    }
	    else {
		track.addLatLng([pt[1], pt[0]]);
	    }
	}
    }

    function showMarkerStep(marker, j) {
	pt = tracklines[marker.track_id].geometry.coordinates[j];
	if (pt !== undefined) {
	    marker.setLatLng([pt[1], pt[0]]);
	}
    }

    function onMapClick(e) {
	console.log("you clicked.", e.latlng);
	var pt = e.latlng;
	url = "http://localhost:8888/drifter?location=" + pt.lat + ',' + pt.lng;
	$.getJSON(url, function(feature) {
	    feature["properties"]["track_id"] = tracklines.length;
	    tracklines.push(feature);
	    var polyline = L.polyline([], {color: 'red', smoothFactor: 0.0});
	    polyline.track_id = tracklines.length - 1;
	    tracks.addLayer(polyline);
	    var marker = L.marker([pt.lat, pt.lng], {icon: duckIcon});
	    marker.track_id = tracklines.length - 1;
	    markers.addLayer(marker);
	});
    }
    var tracks = L.layerGroup([])
    var markers = L.layerGroup([])
    var duckIcon = L.icon({
	iconUrl: 'duck.png',
	iconSize: [16, 16],
    });
    var map = new L.Map('map', {center: latlng, zoom: 7, layers: [basemap, tracks, markers]});
    var tracklines = new Array();
    map.on('click', onMapClick);
    map.addControl(new MyButton({layer: tracks, lines: tracklines, markers: markers}));
    var url = "txla10day.json";

    $.getJSON("domain.json", function(domain) {
	console.log(domain);
	L.geoJson(domain, {style: {fill: false}
			  }).addTo(map);
    })

    animateLines();
  };

MyButton = L.Control.extend({
    // define the "clear tracks" button
    options: {
	position: 'topleft',
	layer: null
    },
    initialize: function (options) {
	this._button = {};
	this.setButton(options);
	this._layer = options.layer;
	this._tracklines = options.lines;
	this._markers = options.markers;
    },
 
    onAdd: function (map, tracks) {
	this._map = map;
	this._tracks = tracks;
	var container = L.DomUtil.create('div', 'leaflet-control-button');
	this._container = container;
	this._update();
	return this._container;
    },
 
    onRemove: function (map) {
    },
 
    clearPoints: function (e) {
	console.log("Clear!");
	this._layer.clearLayers();
	this._tracklines.length = 0;
	this._markers.clearLayers();
    },

    setButton: function (options) {
	var button = {
	    'text': 'Clear tracks', //string
	    'onClick': this.clearPoints, //callback function
	    'hideText': false, //forced bool
	    'maxWidth': 70, //number
	    'doToggle': false,	//bool
	    'toggleStatus': false	//bool
	};
 
	this._button = button;
	this._update();
    },
    getText: function () {
	return this._button.text;
    },
    destroy: function () {
	this._button = {};
	this._update();
    },
    toggle: function (e) {
	if(typeof e === 'boolean'){
	    this._button.toggleStatus = e;
	}
	else{
	    this._button.toggleStatus = !this._button.toggleStatus;
	}
	this._update();
    },
    _update: function () {
	if (!this._map) {
	    return;
	}


	this._container.innerHTML = '';
	this._makeButton(this._button);
    },

    _makeButton: function (button) {
	var newButton = L.DomUtil.create('div', 'leaflet-buttons-control-button', this._container);
	if(button.toggleStatus)
	    L.DomUtil.addClass(newButton,'leaflet-buttons-control-toggleon');
	if(button.text !== ''){
	    
	    L.DomUtil.create('br','',newButton); //there must be a better way
	    
	    var span = L.DomUtil.create('span', 'leaflet-buttons-control-text', newButton);
	    var text = document.createTextNode(button.text); //is there an L.DomUtil for this?
	    span.appendChild(text);
	    if(button.hideText)
		L.DomUtil.addClass(span,'leaflet-buttons-control-text-hide');
	}
	
	L.DomEvent
	    .addListener(newButton, 'click', L.DomEvent.stop)
	    .addListener(newButton, 'click', button.onClick,this)
	    .addListener(newButton, 'click', this._clicked,this);
	L.DomEvent.disableClickPropagation(newButton);
	return newButton;

    },
    _clicked: function () { //'this' refers to button
	if(this._button.doToggle){
	    if(this._button.toggleStatus) {	//currently true, remove class
		L.DomUtil.removeClass(this._container.childNodes[0],'leaflet-buttons-control-toggleon');
	    }
	    else{
		L.DomUtil.addClass(this._container.childNodes[0],'leaflet-buttons-control-toggleon');
	    }
	    this.toggle();
	}
	return;
    }
});
