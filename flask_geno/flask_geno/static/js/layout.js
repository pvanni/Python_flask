window.onload = function() {
  document.getElementById("accountsButton").addEventListener("click", function() {
    var products = document.querySelectorAll(".account");
    for (var i = 0; i < products.length; i++) {
      if (products[i].style.display === "none") {
        products[i].style.display = "block";
      } else {
        products[i].style.display = "none";
      }
    }
  });
  document.getElementById("analysis_runs_Button").addEventListener("click", function() {
    var products = document.querySelectorAll(".runs");
    for (var i = 0; i < products.length; i++) {
      if (products[i].style.display === "none") {
        products[i].style.display = "block";
      } else {
        products[i].style.display = "none";
      }
    }
  });
};

