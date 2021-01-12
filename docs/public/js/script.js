(function (document) {
  var toggle = document.querySelector(".sidebar-toggle");
  var sidebar = document.querySelector("#sidebar");
  var checkbox = document.querySelector("#sidebar-checkbox");

  document.addEventListener(
    "click",
    function (e) {
      var target = e.target;

      if (
        !checkbox.checked ||
        sidebar.contains(target) ||
        target === checkbox ||
        target === toggle
      )
        return;

      checkbox.checked = false;
    },
    false
  );
})(document);

function openTab(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent");
  console.log(tabcontent.length);
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}

//  MODAL WINDOW FOR IMAGES

// create references to the modal...
var modal = document.getElementById("myModal");
// to all images -- note I'm using a class!
var images = document.getElementsByClassName("myImages");
// the image in the modal
var modalImg = document.getElementById("img01");

// Go through all of the images with our custom class
for (var i = 0; i < images.length; i++) {
  var img = images[i];
  // and attach our click listener for this image.
  img.onclick = function (evt) {
    modal.style.display = "block";
    modalImg.src = this.src;
  };
}

// CLOSE MODAL BUTTON
var span = document.getElementsByClassName("close")[0];

span.onclick = function () {
  modal.style.display = "none";
};

// CLOSE MODAL ON ESC PRESS
document.addEventListener("keydown", function (event) {
  if (event.key === "Escape" && modal.style.display === "block") {
    modal.style.display = "none";
  }
});

// CLOSE MODAL WHEN CLICKED ANYWHERE ELSE
window.onclick = function (event) {
  if (event.target == modal) {
    modal.style.display = "none";
  }
};
