window.onload = function () {

    var cloudmadeUrl = 'http://{s}.tile.cloudmade.com/4fd891040a8a4ecb805c388019e23d46/64082/256/{z}/{x}/{y}.png';

    var basemap = new L.TileLayer(cloudmadeUrl, {maxZoom: 18});
    var latlng = new L.LatLng(26.5, -89.51);
    var delay = 20;		// animation delay (larger is slower)
    var Npoints = 100;		// number of points per track
  
    function animateLine(data) {
	// reconstruct a polyline from a GeoJSON object
	var polyline = L.polyline([], {color: 'red'})
	var time = 0;
	for (i in data.geometry.coordinates) {
	    time += delay;
            setTimeout(function(j, pt) {
		return function() {
		    polyline.addLatLng([pt[1], pt[0]]);
		}
            }(i, data.geometry.coordinates[i]), time);
	}
	return polyline;
    }

    function onMapClick(e) {
	console.log("you clicked.", e.latlng);
	var scale = 0.1;
	var pt = e.latlng;
	var polyline = L.polyline([[pt.lat, pt.lng]], {color: 'red'})
	for (var i=0; i<Npoints; i++) {
	    pt.lat += (Math.random() - 0.4) * scale;
	    pt.lng += (Math.random() - 0.4) * scale;
	    polyline.addLatLng([pt.lat, pt.lng]);
	}
	tracks.addLayer(animateLine(polyline.toGeoJSON()));
    }
    var tracks = L.layerGroup([])
    var map = new L.Map('map', {center: latlng, zoom: 6, layers: [basemap, tracks]});
    map.on('click', onMapClick);
    map.addControl(new MyButton({layer: tracks}));

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
