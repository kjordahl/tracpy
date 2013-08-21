window.onload = function () {

    var cloudmadeUrl = 'http://{s}.tile.cloudmade.com/4fd891040a8a4ecb805c388019e23d46/64082/256/{z}/{x}/{y}.png';

    var basemap = new L.TileLayer(cloudmadeUrl, {maxZoom: 18});
    var latlng = new L.LatLng(27, -94);
    var delay = 80;		// animation delay (larger is slower)
    var Npoints = 151;		// number of points per track
    var isRunning = true;

    function showTimeStep(j) {
	tracks.eachLayer(function (layer) {
	    pt = tracklines[layer.track_id].geometry.coordinates[j];
	    if (pt !== undefined) {
		if (j === 0) {
		    layer.setLatLngs([[pt[1], pt[0]]]);
		}
		else {
		    layer.addLatLng([pt[1], pt[0]]);
		}
	    }
	});
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

    function onMapClick(e) {
	console.log("you clicked.", e.latlng);
	var pt = e.latlng;
	url = "http://localhost:8888/drifter?location=" + pt.lat + ',' + pt.lng;
	$.getJSON(url, function(track) {
	    feature = {"type": "Feature",
		       "geometry": track,
		       "properties": {"track_id": tracklines.length}}
	    tracklines.push(feature);
	    var polyline = L.polyline([], {color: 'red', smoothFactor: 0.0});
	    polyline.track_id = tracklines.length - 1;
	    tracks.addLayer(polyline);
	});
    }
    var tracks = L.layerGroup([])
    var map = new L.Map('map', {center: latlng, zoom: 7, layers: [basemap, tracks]});
    var tracklines = new Array();
    map.on('click', onMapClick);
    map.addControl(new MyButton({layer: tracks, lines: tracklines}));
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
